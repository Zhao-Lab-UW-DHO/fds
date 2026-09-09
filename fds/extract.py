from __future__ import annotations
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Tuple
import pysam
from . import indel_realign as _ir
from .variants import Variant
DEFAULT_DEDUP = 'none'
DEFAULT_MIN_MAPQ = 20
DEFAULT_MIN_BQ = 20
DEFAULT_TLEN_MIN = 1
DEFAULT_TLEN_MAX = 500
UMI_MISSING_FATAL_FRAC = 0.5

def call_snv(read, pos0: int, ref_base: str, alt_base: str, min_bq: int, no_bq=None):
    ref_to_q = {r: q for q, r in read.get_aligned_pairs(matches_only=True) if r is not None and q is not None}
    qpos = ref_to_q.get(pos0)
    if qpos is None:
        return None
    seq, quals = (read.query_sequence, read.query_qualities)
    if seq is None or qpos >= len(seq):
        return None
    if quals is None or qpos >= len(quals):
        if no_bq is not None:
            no_bq[0] += 1
    elif quals[qpos] < min_bq:
        return None
    base = seq[qpos]
    if base == alt_base:
        return 'alt'
    if base == ref_base:
        return 'ref'
    return None

def _realign_call(read, hap, min_mapq: int):
    call = _ir.classify_read(read, hap[0], hap[1], hap[2], hap[3], hap[4], min_mapq=min_mapq)
    if call == 'alt':
        return 'alt'
    if call == 'ref':
        return 'ref'
    return None

def _haplotypes(fa, var: Variant):
    try:
        rw, aw, ws, voff, vspan = _ir.build_haplotypes(fa, var.chrom, var.pos - 1, var.ref, var.alt, var.vtype)
    except ValueError:
        return None
    return (_ir._encode(rw), _ir._encode(aw), ws, voff, vspan)

def process_locus(bam, var: Variant, min_mapq: int, min_bq: int, tlen_min: int, tlen_max: int, fa=None, dedup: str=DEFAULT_DEDUP, warn=None, warn_bq=None, indel_len_correction: bool=True, warn_unpaired=None) -> Tuple[List[Tuple[str, int, int]], str]:
    chrom = var.chrom
    pos0 = var.pos - 1
    vt = var.vtype
    is_snv = vt == 'SNP'
    no_bq = [0]
    n_unpaired = [0]
    del_len = len(var.ref) if vt == 'DEL' else 0
    ins_len = len(var.alt) if vt == 'INS' else 0
    hap = None
    if not is_snv:
        if fa is None:
            raise ValueError('a reference FASTA is required to call %s %s:%d %s>%s -- non-SNVs are classified by realignment against REF/ALT haplotypes' % (vt, var.chrom, var.pos, var.ref, var.alt))
        hap = _haplotypes(fa, var)
        if hap is None:
            return ([], 'event_too_large')
    win = max(20, len(var.src_ref), len(var.src_alt))
    frags: Dict[str, dict] = {}
    try:
        it = bam.fetch(chrom, max(0, var.src_pos - win), var.src_pos + win)
    except ValueError:
        return ([], 'contig_absent')
    for read in it:
        if read.is_unmapped or read.is_secondary or read.is_supplementary or (read.cigartuples is None) or (read.mapping_quality < min_mapq):
            continue
        if not read.is_paired:
            n_unpaired[0] += 1
            continue
        if is_snv and (not read.reference_start <= pos0 < read.reference_end):
            continue
        tl = read.template_length
        tlen = abs(tl)
        if tlen < tlen_min or tlen > tlen_max:
            continue
        if is_snv:
            call = call_snv(read, pos0, var.ref, var.alt, min_bq, no_bq)
        else:
            call = _realign_call(read, hap, min_mapq)
        if call is None:
            continue
        if dedup == 'flag' and read.is_duplicate:
            continue
        fs = read.reference_start if tl >= 0 else read.reference_end - tlen
        umi = read.get_tag('RX') if read.has_tag('RX') else ''
        ent = frags.setdefault(read.query_name, {'calls': set(), 'tlen': tlen, 'fs': fs, 'umi': umi})
        ent['calls'].add(call)
        if tl > 0:
            ent['tlen'], ent['fs'] = (tlen, fs)

    def _corrected(allele: str, flen: int, raw: int):
        if allele == 'alt' and indel_len_correction:
            if vt == 'DEL':
                flen -= del_len
            elif vt == 'INS':
                flen += ins_len
            if flen < tlen_min or flen > tlen_max:
                return None
        return (allele, flen, raw)
    out: List[Tuple[str, int]] = []
    seen = set()
    n_no_umi = 0
    n_umi_seen = 0
    for ent in frags.values():
        if dedup == 'umi':
            n_umi_seen += 1
            if not ent['umi']:
                n_no_umi += 1
            key = (ent['fs'], ent['tlen'], ent['umi'])
            if key in seen:
                continue
            seen.add(key)
        calls = ent['calls']
        if calls == {'alt'}:
            allele = 'alt'
        elif calls == {'ref'}:
            allele = 'ref'
        else:
            continue
        c = _corrected(allele, ent['tlen'], ent['tlen'])
        if c is not None:
            out.append(c)
    if n_umi_seen and warn is not None:
        warn(n_no_umi, n_umi_seen)
    if no_bq[0] and warn_bq is not None:
        warn_bq(no_bq[0])
    if warn_unpaired is not None:
        warn_unpaired(n_unpaired[0], len(out))
    if not out and n_unpaired[0]:
        return (out, 'unpaired_only')
    if not out:
        return (out, 'no_coverage')
    return (out, '')

def extract_sample(bam_path: str, variants: Iterable[Variant], min_mapq: int=DEFAULT_MIN_MAPQ, min_bq: int=DEFAULT_MIN_BQ, tlen_min: int=DEFAULT_TLEN_MIN, tlen_max: int=DEFAULT_TLEN_MAX, exclude_fraglens: Iterable[int]=(), reference: str=None, dedup: str=DEFAULT_DEDUP, indel_len_correction: bool=True, progress=None) -> Dict[str, Dict[str, Dict[int, int]]]:
    _no_umi = {'n': 0, 'loci': 0, 'total': 0}

    def _warn(n, total):
        _no_umi['n'] += n
        _no_umi['total'] += total
        if n:
            _no_umi['loci'] += 1
    _unpaired = {'reads': 0, 'out': 0}

    def _warn_unpaired(n, produced):
        _unpaired['reads'] += n
        _unpaired['out'] += produced
    _no_bq = {'n': 0, 'loci': 0}

    def _warn_bq(n):
        _no_bq['n'] += n
        _no_bq['loci'] += 1
    drop = set((int(x) for x in exclude_fraglens))
    counts: Dict[str, Dict[str, Dict[int, int]]] = {}
    fa = pysam.FastaFile(reference) if reference else None
    try:
        _kw = {'reference_filename': reference} if bam_path.lower().endswith('.cram') else {}
        with pysam.AlignmentFile(bam_path, 'rb', **_kw) as bam:
            for i, var in enumerate(variants):
                per = {'alt': defaultdict(int), 'ref': defaultdict(int)}
                pairs, reason = process_locus(bam, var, min_mapq, min_bq, tlen_min, tlen_max, fa, dedup, _warn, _warn_bq, indel_len_correction, _warn_unpaired)
                for allele, flen, raw in pairs:
                    if raw in drop:
                        continue
                    per[allele][flen] += 1
                counts[var.key] = {'alt': dict(per['alt']), 'ref': dict(per['ref']), 'reason': reason}
                if progress is not None and (i + 1) % 500 == 0:
                    progress(i + 1)
    finally:
        if fa is not None:
            fa.close()
    if _unpaired['reads'] and (not _unpaired['out']):
        print(f"WARNING: every read examined was unpaired ({_unpaired['reads']}) and no "
              f'result was produced; this BAM looks single-end. The empty output is not '
              f'absent coverage.', file=sys.stderr)
    if _no_bq['n']:
        print(f"NOTE: {_no_bq['n']} reads across {_no_bq['loci']} loci carry no QUAL; the quality gate was skipped for those.", file=sys.stderr)
    if _no_umi['total']:
        frac = _no_umi['n'] / _no_umi['total']
        if frac >= UMI_MISSING_FATAL_FRAC:
            raise ValueError(f"--dedup umi on {bam_path}: {_no_umi['n']} of {_no_umi['total']} entries ({frac:.0%}) carry no RX tag, so the key degrades to coordinates alone and distinct entries would be merged. This file does not look like a UMI library. Use --dedup flag if its duplicates are marked, or --dedup none if they are already removed.")
        if _no_umi['n']:
            print(f"WARNING: {_no_umi['n']} of {_no_umi['total']} entries ({frac:.1%}, across {_no_umi['loci']} loci) carry no RX tag; the key degraded to coordinates alone for those.", file=sys.stderr)
    return counts
