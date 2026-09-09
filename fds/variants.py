from __future__ import annotations
import csv
import gzip
import io
import os
import re
import stat
from dataclasses import dataclass, field
from itertools import chain
from typing import Dict, List, Optional, Tuple
_ALIASES = {'chrom': ('chrom', 'chr', 'chromosome', 'chrom_hg38', 'contig', 'seqnames', '#chrom'), 'pos': ('pos', 'position', 'pos_hg38', 'start_position'), 'ref': ('ref', 'reference', 'ref_allele', 'reference_allele'), 'alt': ('alt', 'alternate', 'alt_allele', 'tumor_seq_allele2', 'variant_allele')}
_VALID_ALLELE = re.compile('^[ACGTNacgtn]+$')

@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    vtype: str
    src_pos: int
    src_ref: str
    src_alt: str
    filter: str = ''
    meta: Dict[str, str] = field(default_factory=dict)

    @property
    def key(self) -> str:
        return f'{self.chrom}:{self.src_pos}:{self.src_ref}:{self.src_alt}'

class SkipReason:
    SYMBOLIC = 'symbolic_alt'
    COMPLEX = 'complex_allele'
    MALFORMED = 'malformed'
    FILTERED = 'excluded_filter_tag'

_STREAM_CACHE: Dict[str, bytes] = {}

def _is_regular(path: str) -> bool:
    try:
        return stat.S_ISREG(os.stat(path).st_mode)
    except OSError:
        return True

def _raw_bytes(path: str) -> bytes:
    if path not in _STREAM_CACHE:
        with open(path, 'rb') as fh:
            _STREAM_CACHE[path] = fh.read()
    return _STREAM_CACHE[path]

def _open_text(path: str):
    if not _is_regular(path):
        raw = _STREAM_CACHE.pop(path, None)
        if raw is None:
            raw = _raw_bytes(path)
            _STREAM_CACHE.pop(path, None)
        if raw[:2] == b'\x1f\x8b':
            raw = gzip.decompress(raw)
        return io.StringIO(raw.decode('utf-8'), newline='')
    if path.endswith('.gz'):
        return io.TextIOWrapper(gzip.open(path, 'rb'), encoding='utf-8', newline='')
    return open(path, 'r', encoding='utf-8', newline='')

def normalise_alleles(pos: int, ref: str, alt: str) -> Tuple[int, str, str, str]:
    ref = (ref or '').strip()
    alt = (alt or '').strip()
    if alt == '-':
        if not _VALID_ALLELE.match(ref):
            raise ValueError(SkipReason.MALFORMED)
        return (pos, ref.upper(), '-', 'DEL')
    if ref == '-':
        if not _VALID_ALLELE.match(alt):
            raise ValueError(SkipReason.MALFORMED)
        return (pos, '-', alt.upper(), 'INS')
    if alt.startswith('<') or '[' in alt or ']' in alt:
        raise ValueError(SkipReason.SYMBOLIC)
    if not _VALID_ALLELE.match(ref) or not _VALID_ALLELE.match(alt):
        raise ValueError(SkipReason.MALFORMED)
    ref, alt = (ref.upper(), alt.upper())
    if len(ref) == 1 and len(alt) == 1:
        return (pos, ref, alt, 'SNP')
    if len(ref) == len(alt):
        return (pos, ref, alt, 'MNV')
    if len(ref) != len(alt) and ref[0] == alt[0]:
        n_shared = 0
        for a, b in zip(ref, alt):
            if a != b:
                break
            n_shared += 1
        if len(ref) > len(alt):
            if n_shared != len(alt):
                raise ValueError(SkipReason.COMPLEX)
            return (pos + n_shared, ref[n_shared:], '-', 'DEL')
        else:
            if n_shared != len(ref):
                raise ValueError(SkipReason.COMPLEX)
            return (pos + n_shared - 1, '-', alt[n_shared:], 'INS')
    raise ValueError(SkipReason.COMPLEX)

def _make(chrom: str, pos: int, ref: str, alt: str, filt: str, meta: Dict[str, str]) -> Variant:
    npos, nref, nalt, vt = normalise_alleles(pos, ref, alt)
    return Variant(chrom=chrom, pos=npos, ref=nref, alt=nalt, vtype=vt, src_pos=pos, src_ref=ref, src_alt=alt, filter=filt, meta=meta)

def read_vcf(path: str, sample: Optional[str]=None) -> Tuple[List[Variant], Dict[str, int]]:
    variants: List[Variant] = []
    skipped: Dict[str, int] = {}
    samples: List[str] = []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                cols = line.split('\t')
                samples = cols[9:] if len(cols) > 9 else []
                if sample is not None and sample not in samples:
                    raise ValueError(f"--vcf-sample {sample!r} not in {path}; samples are: {(', '.join(samples) if samples else '(none)')}")
                continue
            f = line.split('\t')
            if len(f) < 5:
                skipped[SkipReason.MALFORMED] = skipped.get(SkipReason.MALFORMED, 0) + 1
                continue
            chrom, pos_s, vid, ref, alt_field = (f[0], f[1], f[2], f[3], f[4])
            filt = f[6] if len(f) > 6 else ''
            try:
                pos = int(pos_s)
            except ValueError:
                skipped[SkipReason.MALFORMED] = skipped.get(SkipReason.MALFORMED, 0) + 1
                continue
            gt = None
            if sample is not None and (not samples):
                raise ValueError(f'--vcf-sample {sample!r} was given but {path} has no #CHROM header line')
            if sample is not None and samples:
                idx = 9 + samples.index(sample)
                gt = f[idx] if len(f) > idx else None
            alts = alt_field.split(',')
            for i, alt in enumerate(alts):
                if gt is not None and (not _gt_carries(gt, i + 1)):
                    continue
                try:
                    variants.append(_make(chrom, pos, ref, alt, filt, {'id': vid, 'info': f[7] if len(f) > 7 else ''}))
                except ValueError as e:
                    r = str(e)
                    skipped[r] = skipped.get(r, 0) + 1
    return (variants, skipped)

def _gt_carries(gt_field: str, allele_index: int) -> bool:
    gt = gt_field.split(':')[0]
    if gt in ('.', './.', '.|.'):
        return False
    return str(allele_index) in re.split('[/|]', gt)

def read_table(path: str, colmap: Optional[Dict[str, str]]=None) -> Tuple[List[Variant], Dict[str, int]]:
    variants: List[Variant] = []
    skipped: Dict[str, int] = {}
    with _open_text(path) as fh:
        for head in fh:
            if head.startswith('#') and (not head.startswith('#CHROM')):
                continue
            break
        else:
            return (variants, skipped)
        delim = ',' if head.count(',') > head.count('\t') else '\t'
        rdr = csv.DictReader(chain([head], fh), delimiter=delim)
        fields = [(c or '').strip().lstrip('\ufeff') for c in rdr.fieldnames or []]
        rdr.fieldnames = fields
        lower = {c.lower(): c for c in fields}
        resolved = {}
        for want, aliases in _ALIASES.items():
            if colmap and want in colmap:
                resolved[want] = colmap[want]
                continue
            for a in aliases:
                if a in lower:
                    resolved[want] = lower[a]
                    break
        missing = [w for w in _ALIASES if w not in resolved]
        if missing:
            raise ValueError(f"{path}: cannot find column(s) for {', '.join(missing)}. Header is: {', '.join(fields)}. Use --col chrom=NAME etc. to map them.")
        filt_col = lower.get('filter')
        for row in rdr:
            row = {(k or '').strip(): (v or '').strip() for k, v in row.items()}
            try:
                pos = int(float(row[resolved['pos']]))
            except (KeyError, ValueError):
                skipped[SkipReason.MALFORMED] = skipped.get(SkipReason.MALFORMED, 0) + 1
                continue
            meta = {k: v for k, v in row.items() if k not in resolved.values()}
            try:
                variants.append(_make(row[resolved['chrom']], pos, row[resolved['ref']], row[resolved['alt']], row.get(filt_col, '') if filt_col else '', meta))
            except ValueError as e:
                r = str(e)
                skipped[r] = skipped.get(r, 0) + 1
    return (variants, skipped)

def _is_bcf(path: str) -> bool:
    try:
        if not _is_regular(path):
            raw = _raw_bytes(path)
            if raw[:3] == b'BCF':
                return True
            if raw[:2] == b'\x1f\x8b':
                return gzip.decompress(raw)[:3] == b'BCF'
            return False
        with open(path, 'rb') as fh:
            head = fh.read(4)
        if head[:3] == b'BCF':
            return True
        if head[:2] == b'\x1f\x8b':
            with gzip.open(path, 'rb') as fh:
                return fh.read(3) == b'BCF'
    except OSError:
        return False
    return False

def read_variants(path: str, vcf_sample: Optional[str]=None, colmap: Optional[Dict[str, str]]=None) -> Tuple[List[Variant], Dict[str, int]]:
    base = os.path.basename(path).lower()
    try:
        return _read_variants(path, base, vcf_sample, colmap)
    finally:
        _STREAM_CACHE.pop(path, None)

def _read_variants(path: str, base: str, vcf_sample: Optional[str], colmap: Optional[Dict[str, str]]) -> Tuple[List[Variant], Dict[str, int]]:
    if _is_bcf(path):
        raise ValueError(f'{path} is binary BCF, which this tool cannot read. Convert it first:\n    bcftools view {path} -Oz -o variants.vcf.gz')
    if base.endswith('.vcf') or base.endswith('.vcf.gz'):
        return read_vcf(path, vcf_sample)
    if vcf_sample:
        raise ValueError('--vcf-sample applies to a VCF; the given variant file is a table')
    return read_table(path, colmap)

def apply_filter_tags(variants: List[Variant], exclude: List[str]) -> Tuple[List[Variant], int]:
    if not exclude:
        return (variants, 0)
    ex = {t.strip() for t in exclude if t.strip()}
    kept, dropped = ([], 0)
    for v in variants:
        tags = {t.strip() for t in re.split('[;,]', v.filter or '') if t.strip()}
        if tags & ex:
            dropped += 1
        else:
            kept.append(v)
    return (kept, dropped)
