import numpy as np
MIN_MARGIN = 2
MAX_MISMATCH_FRAC = 0.1
MIN_OVERLAP = 20
FLANK = 150
MAX_EVENT_LEN = 100
_B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def _encode(seq):
    a = np.full(len(seq), 4, dtype=np.int8)
    for b, v in _B.items():
        a[np.frombuffer(seq.encode().upper(), dtype=np.uint8) == ord(b)] = v
    return a

def build_haplotypes(fa, chrom, pos0, ref, alt, vtype, flank=FLANK):
    ev = len(ref) if vtype == 'DEL' else len(alt) if vtype == 'INS' else len(ref)
    if ev > MAX_EVENT_LEN:
        raise ValueError('event of %d bp exceeds MAX_EVENT_LEN=%d: no read can span it, so this method cannot judge it. Use a junction-contig approach.' % (ev, MAX_EVENT_LEN))
    ws = max(0, pos0 - flank)
    we = pos0 + flank + (len(ref) if vtype == 'DEL' else 1)
    ref_win = fa.fetch(chrom, ws, we).upper()
    i = pos0 - ws
    if vtype == 'DEL':
        alt_win = ref_win[:i] + ref_win[i + len(ref):]
    elif vtype == 'INS':
        alt_win = ref_win[:i + 1] + alt.upper() + ref_win[i + 1:]
    else:
        alt_win = ref_win[:i] + alt.upper() + ref_win[i + len(ref):]
    var_span = len(ref) if vtype == 'DEL' else 1 if vtype == 'INS' else len(ref)
    return (ref_win, alt_win, ws, i, var_span)

def _mismatches(read_arr, hap_arr, offset):
    if offset < 0:
        read_arr = read_arr[-offset:]
        offset = 0
    n = min(len(read_arr), len(hap_arr) - offset)
    if n <= 0:
        return (0, 0)
    r = read_arr[:n]
    h = hap_arr[offset:offset + n]
    return (int(np.count_nonzero(r != h)), n)

def classify_read(read, ref_arr, alt_arr, ws, var_off, var_span, min_mapq=20, min_margin=MIN_MARGIN, max_mm_frac=MAX_MISMATCH_FRAC, min_overlap=MIN_OVERLAP):
    if read.mapping_quality < min_mapq:
        return 'low_mapq'
    seq = read.query_sequence
    if not seq or read.cigartuples is None:
        return None
    lead_clip = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
    q_start_ref = read.reference_start - lead_clip
    offset = q_start_ref - ws
    if offset > var_off or offset + len(seq) < var_off + var_span:
        return None
    ra = _encode(seq)
    mm_ref, n_ref = _mismatches(ra, ref_arr, offset)
    mm_alt, n_alt = _mismatches(ra, alt_arr, offset)
    if min(n_ref, n_alt) < min_overlap:
        return None
    if mm_ref <= mm_alt:
        win, lose, n_win, mm_win = ('ref', mm_alt, n_ref, mm_ref)
    else:
        win, lose, n_win, mm_win = ('alt', mm_ref, n_alt, mm_alt)
    if mm_win > max_mm_frac * n_win:
        return 'ambiguous'
    if lose - mm_win < min_margin:
        return 'ambiguous'
    return win

def classify_locus(bam, fa, chrom, pos0, ref, alt, vtype, min_mapq=20, flank=FLANK, win=None, **kw):
    ref_win, alt_win, ws, var_off, var_span = build_haplotypes(fa, chrom, pos0, ref, alt, vtype, flank)
    ref_arr, alt_arr = (_encode(ref_win), _encode(alt_win))
    span = len(ref) if vtype == 'DEL' else 1
    w = win if win is not None else max(20, span + 5)
    out = []
    for read in bam.fetch(chrom, max(0, pos0 - w), pos0 + w + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or (read.cigartuples is None):
            continue
        call = classify_read(read, ref_arr, alt_arr, ws, var_off, var_span, min_mapq=min_mapq, **kw)
        if call is not None:
            out.append((read, call))
    return out
