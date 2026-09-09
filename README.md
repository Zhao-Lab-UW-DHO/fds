# fds

Per-variant scoring from aligned reads.

Takes a list of variants that have already been identified, reads the supporting and
non-supporting observations for each one out of an indexed BAM or CRAM, and applies a JSON
profile.

## Requirements

Python ≥ 3.8, `pysam`, `numpy`.

## Use

There is no install step: run it from the repository root (or put that directory on
`PYTHONPATH`), as a module.

```bash
# a cohort: a TSV with sample, bam, variants columns
python3 -m fds.cli score --manifest cohort.tsv --profile profile.json \
        --reference ref.fa --out scores.tsv.gz

# one sample
python3 -m fds.cli score --bam s.bam --variants s.vcf --sample S1 \
        --profile profile.json --reference ref.fa --out scores.tsv
```

`--reference` is an indexed FASTA matching the alignment's build. It is required when the
variant list contains a non-SNV, and always for CRAM input.

## The profile

`fds` applies a locked model supplied as a JSON file (`--profile`). **No profile is distributed
with this repository.** The profile is published as a Supplementary Table with the associated
manuscript; download and pass the json file and use its path with `--profile`.

A profile is validated when it is loaded and refused if its own parameters disagree with each
other, so a file damaged or truncated in transit fails loudly instead of scoring.

## Inputs

Variants may be supplied as VCF, TSV, CSV or MAF. Four fields are required — chromosome,
position, reference allele, alternate allele — under any of the usual column spellings; use
`--col` to map a non-standard header. Both anchor-base and MAF indel encodings are accepted.


## Output

One row per variant. Nothing that can be read is dropped: a variant that cannot be evaluated is
flagged rather than omitted, and `reason` says why — `no_coverage`, `contig_absent`,
`event_too_large` or `unpaired_only`, empty when observations were found. A locus with no
observations at all scores `NA` in every score column rather than a numeric zero.

Two things do leave the table, both by request and both reported on stderr: records the four-field
contract cannot represent (symbolic alleles, breakends, complex substitutions) are skipped and
tallied, and `--exclude-filter-tags` removes rows when you ask it to.

`--keep-source-columns` carries every other column of the variant file through, except a VCF's
`INFO` and any column whose name collides with a computed one — the run names what it dropped.

⚠️ The calibrated columns **rank** variants. Their absolute values are not probabilities for any
particular cohort, and should not be read as such.

## Tests

```bash
python3 -m pytest -q
```
