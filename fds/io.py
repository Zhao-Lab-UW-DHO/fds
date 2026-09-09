from __future__ import annotations
import csv
import gzip
import io
import os
from typing import Dict, List, Optional

def _fmt(v, digits: int=6) -> str:
    if v is None:
        return 'NA'
    if isinstance(v, bool):
        return 'TRUE' if v else 'FALSE'
    if isinstance(v, float):
        if v != v:
            return 'NA'
        return f'{v:.{digits}g}'
    return str(v)

class TableWriter:

    def __init__(self, path: str, columns: List[str], digits: int=6):
        self.columns = list(columns)
        self.digits = digits
        self.path = path
        if path.endswith('.gz'):
            self._raw = gzip.open(path, 'wb')
            self._fh = io.TextIOWrapper(self._raw, encoding='utf-8', newline='')
        else:
            self._raw = None
            self._fh = open(path, 'w', encoding='utf-8', newline='')
        self._w = csv.writer(self._fh, delimiter='\t', lineterminator='\n')
        self._w.writerow(self.columns)

    def write(self, row: Dict[str, object]) -> None:
        self._w.writerow([_fmt(row.get(c), self.digits) for c in self.columns])

    def close(self) -> None:
        self._fh.close()
        if self._raw is not None:
            self._raw.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

def ensure_case_sensitive(path: str) -> Optional[str]:
    d = path if os.path.isdir(path) else os.path.dirname(path) or '.'
    probe_lower = os.path.join(d, '.fds_case_probe_a')
    probe_upper = os.path.join(d, '.fds_case_probe_A')
    try:
        with open(probe_lower, 'w') as fh:
            fh.write('x')
        collides = os.path.exists(probe_upper)
        os.unlink(probe_lower)
        if collides:
            return f'{d} is on a CASE-INSENSITIVE filesystem. Samples whose ids differ only in case (e.g. 733-A vs 733-a) will overwrite each other.'
    except OSError:
        return None
    return None
