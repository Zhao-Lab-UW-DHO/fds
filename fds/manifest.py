from __future__ import annotations
import csv
import os
from dataclasses import dataclass
from typing import Dict, List

@dataclass
class Job:
    sample: str
    bam: str
    variants: str

def read_manifest(path: str) -> List[Job]:
    jobs: List[Job] = []
    seen: Dict[str, int] = {}
    with open(path, 'r', encoding='utf-8', newline='') as fh:
        head = fh.readline()
        if not head:
            raise ValueError(f'{path}: empty manifest')
        delim = ',' if head.count(',') > head.count('\t') else '\t'
        fh.seek(0)
        rdr = csv.DictReader(fh, delimiter=delim)
        cols = {(c or '').strip().lower(): (c or '').strip() for c in rdr.fieldnames or []}
        for req in ('sample', 'bam', 'variants'):
            if req not in cols:
                raise ValueError(f"{path}: manifest needs a {req!r} column; header is {', '.join(cols.values()) or '(empty)'}")
        for i, row in enumerate(rdr, start=2):
            row = {k: (v or '').strip() for k, v in row.items()}
            sample = row[cols['sample']]
            bam = row[cols['bam']]
            variants = row[cols['variants']]
            missing = [n for n, v in (('sample', sample), ('bam', bam), ('variants', variants)) if not v]
            if missing:
                raise ValueError(f"{path} line {i}: {', '.join(missing)} must not be empty")
            if sample in seen:
                raise ValueError(f'{path} line {i}: sample {sample!r} already appears on line {seen[sample]}. This tool takes ONE BAM and ONE variant list per sample -- combine them upstream (samtools merge / bcftools concat) and give it the merged file.')
            seen[sample] = i
            jobs.append(Job(sample=sample, bam=bam, variants=variants))
    if not jobs:
        raise ValueError(f'{path}: manifest has no rows')
    return jobs

def preflight(jobs: List[Job], reference: str=None) -> List[str]:
    problems: List[str] = []
    if reference:
        if not os.path.exists(reference):
            problems.append(f'reference FASTA not found: {reference}')
        elif not os.path.exists(reference + '.fai'):
            problems.append(f'reference FASTA index not found beside {reference} (run: samtools faidx {reference})')
    for j in jobs:
        if not os.path.exists(j.bam):
            problems.append(f'{j.sample}: BAM not found: {j.bam}')
        elif not any((os.path.exists(p) for p in (j.bam + '.bai', os.path.splitext(j.bam)[0] + '.bai', j.bam + '.csi', j.bam + '.crai', os.path.splitext(j.bam)[0] + '.crai'))):
            problems.append(f'{j.sample}: alignment index not found beside {j.bam} (.bai/.csi for BAM, .crai for CRAM)')
        elif j.bam.lower().endswith('.cram') and (not reference):
            problems.append(f'{j.sample}: {j.bam} is CRAM, which cannot be decoded without --reference (the FASTA it was compressed against)')
        if not os.path.exists(j.variants):
            problems.append(f'{j.sample}: variants file not found: {j.variants}')
    return problems
