from __future__ import annotations
import math
from typing import Dict, Optional
import numpy as np
from .profile import Profile

def _bin_index(flen: int, lo: int, hi: int) -> int:
    if flen < lo or flen > hi:
        return -1
    return flen - lo

def shares(counts: Dict[int, int], lo: int, hi: int) -> np.ndarray:
    n_bins = hi - lo + 1
    v = np.zeros(n_bins, dtype=float)
    total = 0
    for flen, n in counts.items():
        b = _bin_index(int(flen), lo, hi)
        if b >= 0:
            v[b] += n
            total += n
    if total > 0:
        v /= total
    else:
        v[:] = 0.0
    return v

def movsum(v: np.ndarray, width: int) -> np.ndarray:
    if v.size < width:
        raise ValueError(f'cannot take a {width}-wide moving sum of {v.size} bins')
    cs = np.concatenate(([0.0], np.cumsum(v)))
    return cs[width:] - cs[:-width]

def design_row(alt_counts: Dict[int, int], ref_counts: Dict[int, int], prof: Profile) -> np.ndarray:
    lo, hi = (prof.score_tlen_min, prof.score_tlen_max)
    p = movsum(shares(alt_counts, lo, hi), prof.window_width)
    w = movsum(shares(ref_counts, lo, hi), prof.window_width)
    return p - w

def vaf_of(alt_counts: Dict[int, int], ref_counts: Dict[int, int], prof: Profile) -> Optional[float]:
    lo, hi = (prof.vaf_frag_min, prof.vaf_frag_max)
    a = sum((n for f, n in alt_counts.items() if lo <= int(f) <= hi))
    r = sum((n for f, n in ref_counts.items() if lo <= int(f) <= hi))
    if a + r <= 0:
        return None
    return a / (a + r)

def n_in_window(counts: Dict[int, int], lo: int, hi: int) -> int:
    return sum((n for f, n in counts.items() if lo <= int(f) <= hi))

def calibrate(score: Optional[float], intercept: float, slope: float) -> Optional[float]:
    if score is None or not math.isfinite(score):
        return None
    z = intercept + slope * score
    if z >= 0:
        return 1.0 / (1.0 + math.exp(-z))
    e = math.exp(z)
    return e / (1.0 + e)

class Scorer:

    def __init__(self, prof: Profile):
        self.prof = prof

    def score_variant(self, alt_counts: Dict[int, int], ref_counts: Dict[int, int]) -> Dict[str, object]:
        p = self.prof
        x = design_row(alt_counts, ref_counts, p)
        vaf = vaf_of(alt_counts, ref_counts, p)
        n_alt = n_in_window(alt_counts, p.score_tlen_min, p.score_tlen_max)
        n_ref = n_in_window(ref_counts, p.score_tlen_min, p.score_tlen_max)
        out: Dict[str, object] = {'n_alt': n_alt, 'n_ref': n_ref, 'n_alt_vaf_window': n_in_window(alt_counts, p.vaf_frag_min, p.vaf_frag_max), 'n_ref_vaf_window': n_in_window(ref_counts, p.vaf_frag_min, p.vaf_frag_max), 'vaf': vaf, 'gate_pass': n_alt >= p.gate_min_alt_frags, 'has_ref': n_ref > 0}
        if n_alt == 0 and n_ref == 0:
            for arm in p.arms.values():
                out[arm.out_score] = None
                out[arm.out_prob] = None
            return out
        for arm_name, arm in p.arms.items():
            s = float(np.dot(x, arm.beta))
            if arm.uses_vaf:
                if vaf is None:
                    out[arm.out_score] = None
                    out[arm.out_prob] = None
                    continue
                s = s + arm.beta_vaf * vaf
            out[arm.out_score] = s
            out[arm.out_prob] = calibrate(s, arm.cal_intercept, arm.cal_slope)
        return out
