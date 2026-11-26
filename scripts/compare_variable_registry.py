"""Compare every shared CSV/Excel variable field by proposed_name."""
from __future__ import annotations

import csv
from decimal import Decimal, InvalidOperation
from pathlib import Path

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "docs" / "variable_registry.csv"
XLSX_PATH = ROOT / "docs" / "HTGR_전체시스템_통합변수표_SG_RCCS_보완반영.xlsx"
EXPECTED_COUNT = 460
FIELDS = (
    ("component", 0), ("proposed_name", 1), ("variable", 2),
    ("kind_tentative", 3), ("unit", 4), ("shape", 5),
    ("stage", 6), ("reference_value", 7), ("source", 8),
    ("status", 9), ("notes", 10), ("original_notation", 11),
)


def _text(value):
    return "" if value is None else str(value).strip()


def _equal(field, csv_value, excel_value):
    if field != "reference_value":
        return _text(csv_value) == _text(excel_value)
    if _text(csv_value) == _text(excel_value):
        return True
    def normalized(value):
        if isinstance(value, bool):
            return Decimal(int(value))
        spelling = _text(value)
        if spelling.casefold() in ("true", "false"):
            return Decimal(int(spelling.casefold() == "true"))
        try:
            return Decimal(spelling)
        except (InvalidOperation, ValueError, TypeError):
            return spelling
    return normalized(csv_value) == normalized(excel_value)


def compare(csv_path=CSV_PATH, xlsx_path=XLSX_PATH):
    """Return all mismatches as (name, field, csv_value, excel_value)."""
    with Path(csv_path).open(encoding="utf-8-sig", newline="") as handle:
        csv_rows = list(csv.DictReader(handle))
    workbook = load_workbook(xlsx_path, read_only=True, data_only=True)
    try:
        excel_rows = list(workbook["전체 변수표"].iter_rows(min_row=2, values_only=True))
    finally:
        workbook.close()
    csv_by_name = {_text(row["proposed_name"]): row for row in csv_rows}
    excel_by_name = {_text(row[1]): row for row in excel_rows}
    mismatches = []
    if len(csv_by_name) != len(csv_rows):
        mismatches.append(("<CSV>", "duplicate proposed_name", len(csv_rows), len(csv_by_name)))
    if len(excel_by_name) != len(excel_rows):
        mismatches.append(("<Excel>", "duplicate proposed_name", len(excel_rows), len(excel_by_name)))
    for name in sorted(csv_by_name.keys() | excel_by_name.keys()):
        if name not in csv_by_name:
            mismatches.append((name, "missing CSV row", None, "present"))
            continue
        if name not in excel_by_name:
            mismatches.append((name, "missing Excel row", "present", None))
            continue
        csv_row, excel_row = csv_by_name[name], excel_by_name[name]
        for field, index in FIELDS:
            if not _equal(field, csv_row[field], excel_row[index]):
                mismatches.append((name, field, csv_row[field], excel_row[index]))
    return len(csv_rows), len(excel_rows), mismatches


if __name__ == "__main__":
    count_csv, count_excel, differences = compare()
    print(f"CSV={count_csv}, Excel={count_excel}, semantic mismatches={len(differences)}")
    for name, field, left, right in differences:
        print(f"{name} / {field}: CSV={left!r}; Excel={right!r}")
    raise SystemExit(bool(differences or count_csv != EXPECTED_COUNT or count_excel != EXPECTED_COUNT))
