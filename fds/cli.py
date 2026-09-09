from __future__ import annotations
import argparse
import json
import os
import re
import sys
from typing import Dict, List, Optional
from . import __version__
from .extract import DEFAULT_DEDUP, DEFAULT_MIN_BQ, DEFAULT_MIN_MAPQ, extract_sample
from .io import TableWriter, ensure_case_sensitive
from .manifest import Job, preflight, read_manifest
from .profile import load_profile
from .score import Scorer
from .variants import apply_filter_tags, read_variants
BASE_COLUMNS = ['sample', 'chrom', 'pos', 'ref', 'alt', 'variant_type', 'n_alt', 'n_ref', 'n_alt_vaf_window', 'n_ref_vaf_window', 'vaf', 'gate_pass', 'has_ref', 'reason', 'filter']

def _parse_colmap(pairs: List[str]) -> Dict[str, str]:
    out = {}
    for p in pairs or []:
        if '=' not in p:
            raise SystemExit(f'--col expects name=COLUMN, got {p!r}')
        k, v = p.split('=', 1)
        k = k.strip().lower()
        if k not in ('chrom', 'pos', 'ref', 'alt'):
            raise SystemExit(f'--col name must be one of chrom/pos/ref/alt, got {k!r}')
        out[k] = v.strip()
    return out

def _load_profile(path: str):
    try:
        return load_profile(path)
    except FileNotFoundError:
        raise SystemExit(f'profile not found: {path}')
    except IsADirectoryError:
        raise SystemExit(f'--profile expects a JSON file, but {path} is a directory')
    except json.JSONDecodeError as e:
        raise SystemExit(f'{path}: not valid JSON ({e})')
    except ValueError as e:
        raise SystemExit(str(e))
    except KeyError as e:
        raise SystemExit(f'{path}: profile is missing {e} -- it is not a complete fds profile')
    except (TypeError, IndexError) as e:
        raise SystemExit(f'{path}: profile is malformed ({type(e).__name__}: {e})')
    except OSError as e:
        raise SystemExit(f'{path}: cannot be read ({e})')


def cmd_score(args) -> int:
    prof = _load_profile(args.profile)
    if args.manifest:
        if args.bam or args.variants or args.sample:
            raise SystemExit('--manifest is exclusive with --bam / --variants / --sample')
        try:
            jobs = read_manifest(args.manifest)
        except ValueError as e:
            raise SystemExit(str(e))
    else:
        if not (args.bam and args.variants):
            raise SystemExit('give --manifest, or both --bam and --variants')
        sample = args.sample or re.sub('\\.(bam|cram|sam)(\\.gz)?$', '', os.path.basename(args.bam), flags=re.I)
        jobs = [Job(sample=sample, bam=args.bam, variants=args.variants)]
    problems = preflight(jobs, args.reference)
    if problems:
        for p in problems:
            print(f'  {p}', file=sys.stderr)
        raise SystemExit(f'{len(problems)} missing input(s); nothing was run.')
    out_dir = os.path.dirname(os.path.abspath(args.out))
    if not os.path.isdir(out_dir):
        raise SystemExit(f'output directory does not exist: {out_dir}')
    if not os.access(out_dir, os.W_OK):
        raise SystemExit(f'output directory is not writable: {out_dir}')
    warn = ensure_case_sensitive(args.out)
    if warn:
        print(f'WARNING: {warn}', file=sys.stderr)
    colmap = _parse_colmap(args.col)
    exclude_tags = [t for t in (args.exclude_filter_tags or '').split(',') if t.strip()]
    exclude_fraglens = []
    for x in (args.exclude_fraglens or '').split(','):
        if not x.strip():
            continue
        try:
            exclude_fraglens.append(int(x))
        except ValueError:
            raise SystemExit(f'--exclude-fraglens takes whole numbers separated by commas; {x.strip()!r} is not one.')
    scorer = Scorer(prof)
    arm_cols: List[str] = []
    for a in prof.arms.values():
        arm_cols += [a.out_score, a.out_prob]
    passthrough: List[str] = []
    shadowed: List[str] = []
    per_job = []
    for j in jobs:
        try:
            vs, skipped = read_variants(j.variants, args.vcf_sample, colmap)
        except ValueError as e:
            raise SystemExit(f'[{j.sample}] {e}')
        vs, n_tag_dropped = apply_filter_tags(vs, exclude_tags)
        per_job.append((j, vs, skipped, n_tag_dropped))
        for v in vs:
            for k in v.meta:
                _reserved = {c.lower() for c in BASE_COLUMNS} | {c.lower() for c in arm_cols}
                if k in passthrough or k in ('info',) or k in shadowed:
                    continue
                if k.lower() in _reserved:
                    shadowed.append(k)
                    continue
                passthrough.append(k)
    if not args.reference:
        for j, vs, _s, _t in per_job:
            bad = next((v for v in vs if v.vtype != 'SNP'), None)
            if bad is not None:
                raise SystemExit(f'[{j.sample}] --reference is required: this variant list contains non-SNVs (e.g. {bad.chrom}:{bad.pos} {bad.ref}>{bad.alt}, {bad.vtype}), which are called by realignment against REF/ALT haplotypes built from the genome. Pass the FASTA the BAMs were aligned to.')
    if args.keep_source_columns and shadowed:
        print(f"WARNING: {len(shadowed)} source column(s) were NOT carried through because this "
              f"tool computes a column of the same name: {', '.join(sorted(shadowed))}. "
              f"The values in the output are this tool's, not the variant file's. "
              f'Rename them upstream if you need both.', file=sys.stderr)
    if not args.keep_source_columns:
        passthrough = []
    columns = BASE_COLUMNS + arm_cols + passthrough
    n_written = 0
    with TableWriter(args.out, columns) as w:
        for j, variants, skipped, n_tag_dropped in per_job:
            msg = f'[{j.sample}] {len(variants)} variants'
            if skipped:
                msg += ' | skipped ' + ', '.join((f'{k}={v}' for k, v in sorted(skipped.items())))
            if n_tag_dropped:
                msg += f' | filter-tag dropped {n_tag_dropped}'
            print(msg, file=sys.stderr)
            if not variants:
                continue
            try:
                counts = extract_sample(j.bam, variants, min_mapq=args.min_mapq, min_bq=args.min_base_qual, tlen_min=prof.vaf_frag_min, tlen_max=prof.vaf_frag_max, exclude_fraglens=exclude_fraglens, reference=args.reference, dedup=args.dedup, indel_len_correction=prof.indel_len_correction, progress=(lambda n, s=j.sample: print(f'  [{s}] {n} loci', file=sys.stderr)) if args.verbose else None)
            except ValueError as e:
                raise SystemExit(f'[{j.sample}] {e}')
            n_gated = 0
            for v in variants:
                c = counts.get(v.key, {'alt': {}, 'ref': {}, 'reason': ''})
                res = scorer.score_variant(c['alt'], c['ref'])
                row = {'sample': j.sample, 'chrom': v.chrom, 'pos': v.src_pos, 'ref': v.src_ref, 'alt': v.src_alt, 'variant_type': v.vtype, 'reason': c.get('reason', ''), 'filter': v.filter}
                row.update(res)
                if args.keep_source_columns:
                    row.update({k: val for k, val in v.meta.items() if k in passthrough})
                w.write(row)
                n_written += 1
                if res['gate_pass']:
                    n_gated += 1
            print(f'  [{j.sample}] {n_gated}/{len(variants)} passed the >={prof.gate_min_alt_frags} gate', file=sys.stderr)
    print(f'wrote {args.out}: {n_written} rows over {len(jobs)} sample(s)', file=sys.stderr)
    return 0

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog='fds', description='Per-variant scoring from aligned reads.')
    ap.add_argument('--version', action='version', version=f'fds {__version__}')
    sub = ap.add_subparsers(dest='cmd', required=True)
    sc = sub.add_parser('score', help='score variants')
    sc.add_argument('--profile', required=True, help='profile JSON')
    sc.add_argument('--out', required=True, help='output TSV (.gz allowed)')
    sc.add_argument('--manifest', help='TSV with sample, bam, variants -- all three required, one row per sample. A repeated sample id is an error: merge upstream instead.')
    sc.add_argument('--bam', help='single-sample mode: indexed BAM')
    sc.add_argument('--variants', help='single-sample mode: variant file')
    sc.add_argument('--sample', help='single-sample mode: sample id')
    sc.add_argument('--vcf-sample', help="OPT-IN: subset a multi-sample VCF to records this sample's GT carries. Default is to read the VCF as a plain locus list.")
    sc.add_argument('--col', action='append', default=[], help='map a variant-table column, e.g. --col chrom=Chromosome (repeatable; chrom/pos/ref/alt only)')
    sc.add_argument('--exclude-filter-tags', default='', help='comma-separated FILTER tags to drop; empty by default')
    sc.add_argument('--exclude-fraglens', default='', help='comma-separated values to exclude')
    sc.add_argument('--reference', default=None, help='indexed reference FASTA')
    sc.add_argument('--dedup', choices=('none', 'umi', 'flag'), default=DEFAULT_DEDUP, help='duplicate handling for this BAM')
    sc.add_argument('--min-mapq', type=int, default=DEFAULT_MIN_MAPQ)
    sc.add_argument('--min-base-qual', type=int, default=DEFAULT_MIN_BQ)
    sc.add_argument('--keep-source-columns', action='store_true', help='carry every other column of the variant file into the output')
    sc.add_argument('--verbose', action='store_true')
    sc.set_defaults(func=cmd_score)
    return ap

def main(argv: Optional[List[str]]=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)
if __name__ == '__main__':
    sys.exit(main())
