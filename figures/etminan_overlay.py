"""Load and plot Etminan phi-alpha reference curves from validation CSV."""
from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ETMINAN_CSV = ROOT / "validation" / "etminan_phialpha_reference.csv"


def read_etminan_csv(path: Path) -> dict[str, list[float]] | None:
    if not path.exists():
        return None
    with path.open(newline="") as f:
        lines = [
            line for line in f
            if line.strip() and not line.lstrip().startswith("#")
        ]
    rows = list(csv.DictReader(lines))
    if not rows or "k_MeV" not in rows[0]:
        return None
    out: dict[str, list[float]] = {h: [] for h in rows[0].keys()}
    for row in rows:
        for key, value in row.items():
            if value is None or str(value).strip() == "":
                out[key].append(float("nan"))
            else:
                out[key].append(float(value))
    return out


def _column_for_R(data: dict[str, list[float]], R: float) -> str | None:
    rs = f"R{R:.2f}"
    for name in (f"C_{rs}", f"C_{R:g}", f"C_R{R:.2f}", f"C_R{R:g}"):
        if name in data and not all(math.isnan(v) for v in data[name]):
            return name
    for name in ("C_center", "C", "C_KP", "C_LL"):
        if name in data and not all(math.isnan(v) for v in data[name]):
            return name
    for key in data:
        if key == "k_MeV":
            continue
        if not all(math.isnan(v) for v in data[key]):
            return key
    return None


def etminan_column_for_R(data: dict[str, list[float]] | None, R: float) -> str | None:
    if not data:
        return None
    return _column_for_R(data, R)


def plot_etminan_reference(ax, data: dict[str, list[float]] | None, R: float, *, label: str = "Etminan published reference") -> bool:
    col = etminan_column_for_R(data, R)
    if not data or not col:
        return False
    ax.plot(data["k_MeV"], data[col], color="0.25", lw=2.0, label=label)
    lower = f"{col}_lower" if f"{col}_lower" in data else "C_lower"
    upper = f"{col}_upper" if f"{col}_upper" in data else "C_upper"
    if lower in data and upper in data:
        ax.fill_between(
            data["k_MeV"],
            data[lower],
            data[upper],
            color="0.25",
            alpha=0.15,
            lw=0,
        )
    return True


def etminan_caption_note(data: dict[str, list[float]] | None) -> str:
    if not data:
        return "Etminan published curve overlay: pending (add validation/etminan_phialpha_reference.csv)."
    cols = {R: etminan_column_for_R(data, R) for R in (1.0, 3.0, 5.0)}
    if len(set(cols.values())) > 1:
        return "Etminan reference CSV loaded with per-R columns; qualitative comparison only."
    return "Etminan reference CSV loaded; qualitative comparison only."
