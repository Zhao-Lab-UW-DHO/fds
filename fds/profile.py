from __future__ import annotations
import json
from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np
_ARM_COLUMNS = {'FDS': ('fds', 'fds_p'), 'FDS_VAF': ('fds_vaf', 'fds_vaf_p')}

def _scalar_or_none(v):
    if v is None:
        return None
    if isinstance(v, (dict, list, tuple)) and len(v) == 0:
        return None
    if isinstance(v, (list, tuple)) and len(v) == 1:
        return float(v[0])
    return float(v)

@dataclass
class Arm:
    name: str
    beta: np.ndarray
    beta_vaf: Optional[float]
    uses_vaf: bool
    cal_intercept: float
    cal_slope: float
    out_score: str
    out_prob: str
    train_pooled_auc: Optional[float] = None
    cal_note: Optional[str] = None

@dataclass
class Profile:
    model: str
    source_script: str
    exported: str
    gate_min_alt_frags: int
    score_tlen_min: int
    score_tlen_max: int
    window_width: int
    vaf_frag_min: int
    vaf_frag_max: int
    seed: Optional[int]
    indel_len_correction: bool
    n_windows: int
    arms: Dict[str, Arm]
    arm_document: str = '?'

    @property
    def n_bins(self) -> int:
        return self.score_tlen_max - self.score_tlen_min + 1


def load_profile(path: str) -> Profile:
    with open(path, 'r', encoding='utf-8') as fh:
        d = json.load(fh)
    sv = d.get('schema_version')
    if sv is not None and int(sv) < 2:
        raise ValueError(f'{path}: unsupported schema_version {sv}')
    k = d.get('knobs')
    if not k:
        raise ValueError(f'{path}: no `knobs` block -- this is not an fds profile')
    for req in ('gate_min_alt_frags', 'score_tlen_min', 'score_tlen_max', 'window_width', 'vaf_frag_min', 'vaf_frag_max'):
        if req not in k:
            raise ValueError(f'{path}: knobs is missing {req!r}')
    if k['score_tlen_min'] < k['vaf_frag_min'] or k['score_tlen_max'] > k['vaf_frag_max']:
        raise ValueError(f"{path}: knobs are inconsistent -- [score_tlen_min, score_tlen_max] must lie inside [vaf_frag_min, vaf_frag_max]; got [{k['score_tlen_min']}, {k['score_tlen_max']}] and [{k['vaf_frag_min']}, {k['vaf_frag_max']}].")
    arm_document = d.get('arm_document', '') or ''
    windows = d.get('windows') or {}
    lo = windows.get('lo') or []
    hi = windows.get('hi') or []
    expected = k['score_tlen_max'] - k['score_tlen_min'] + 1 - k['window_width'] + 1
    if expected < 1:
        raise ValueError(f"{path}: knobs are inconsistent -- window_width {k['window_width']} exceeds the span [{k['score_tlen_min']}, {k['score_tlen_max']}], so no window fits.")
    n_windows = len(lo) or expected
    geom = k.get('window_geometry')
    if geom is not None and geom != 'overlap':
        raise ValueError(f'{path}: unsupported window_geometry {geom!r}')
    if lo and n_windows != expected:
        raise ValueError(f'{path}: profile declares {n_windows} windows but its own knobs imply {expected}.')
    arms: Dict[str, Arm] = {}
    for name, a in (d.get('arms') or {}).items():
        beta = np.asarray(a['beta'], dtype=float)
        if lo:
            if lo[0] != k['score_tlen_min']:
                raise ValueError(f'{path}: windows start at {lo[0]}, not {k["score_tlen_min"]}')
            if any(lo[i + 1] - lo[i] != 1 for i in range(len(lo) - 1)):
                raise ValueError(f'{path}: window starts are not stride-1 contiguous')
            if hi and any(h - l != k['window_width'] for l, h in zip(lo, hi)):
                raise ValueError(f'{path}: hi - lo is not window_width {k["window_width"]}')
        if beta.size != n_windows:
            raise ValueError(f'{path}: arm {name} has {beta.size} coefficients but the profile declares {n_windows} windows')
        cal = a.get('cal_row_score')
        if not cal:
            raise ValueError(f'{path}: arm {name} is missing cal_row_score')
        uses_vaf = bool(a.get('uses_vaf', False))
        bvaf = _scalar_or_none(a.get('beta_vaf'))
        if uses_vaf and bvaf is None:
            raise ValueError(f'{path}: arm {name} declares uses_vaf but has no beta_vaf')
        cols = _ARM_COLUMNS.get(name, (name.lower(), f'{name.lower()}_p'))
        arms[name] = Arm(name=name, beta=beta, beta_vaf=bvaf, uses_vaf=uses_vaf, cal_intercept=float(cal['intercept']), cal_slope=float(cal['score']), out_score=cols[0], out_prob=cols[1], train_pooled_auc=a.get('train_pooled_auc'), cal_note=cal.get('note'))
    if not arms:
        raise ValueError(f'{path}: no arms in profile')
    return Profile(model=d.get('model', '?'), source_script=d.get('source_script', '?'), arm_document=arm_document or '?', exported=d.get('exported', '?'), gate_min_alt_frags=int(k['gate_min_alt_frags']), score_tlen_min=int(k['score_tlen_min']), score_tlen_max=int(k['score_tlen_max']), window_width=int(k['window_width']), vaf_frag_min=int(k['vaf_frag_min']), vaf_frag_max=int(k['vaf_frag_max']), seed=k.get('seed'), indel_len_correction=bool(k.get('indel_len_correction', True)), n_windows=n_windows or expected, arms=arms)
