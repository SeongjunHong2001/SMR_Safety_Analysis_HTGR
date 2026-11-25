import csv

from scripts.compare_variable_registry import CSV_PATH, XLSX_PATH, _equal, compare


def test_all_372_shared_variable_fields_match_csv_and_excel():
    csv_count, excel_count, differences = compare()
    assert csv_count == excel_count == 372
    assert differences == []


def test_compare_detects_semantic_drift_and_normalizes_boolean(tmp_path):
    assert _equal("reference_value", "TRUE", True)
    assert _equal("reference_value", "0", False)
    assert _equal("reference_value", "1e-3", 0.001)
    with CSV_PATH.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["notes"] = "changed description"
    rows[0]["status"] = "changed status"
    changed = tmp_path / "registry.csv"
    with changed.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    _, _, differences = compare(changed, XLSX_PATH)
    assert {(name, field) for name, field, *_ in differences} == {
        (rows[0]["proposed_name"], "notes"),
        (rows[0]["proposed_name"], "status"),
    }
