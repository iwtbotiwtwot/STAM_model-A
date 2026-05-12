"""Quick read of Sean's original SU spreadsheets to inspect derivation."""

import sys
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import openpyxl

files = [
    r"C:\Users\drwho\OneDrive\Documents\SU625\STAM.xlsx",
    r"C:\Users\drwho\OneDrive\Documents\SU625\SU625.xlsx",
]

for fpath in files:
    print("=" * 78)
    print(f"FILE: {fpath}")
    print("=" * 78)
    try:
        wb = openpyxl.load_workbook(fpath, data_only=False)
        for sheet_name in wb.sheetnames:
            ws = wb[sheet_name]
            print(f"\n--- Sheet: {sheet_name} ---")
            print(f"Dimensions: {ws.dimensions}, rows={ws.max_row}, cols={ws.max_column}")

            # Print first 60 rows, all populated columns
            max_print_rows = min(ws.max_row, 60)
            max_print_cols = min(ws.max_column, 12)
            for row_idx in range(1, max_print_rows + 1):
                row_vals = []
                for col_idx in range(1, max_print_cols + 1):
                    cell = ws.cell(row=row_idx, column=col_idx)
                    v = cell.value
                    if v is None:
                        row_vals.append("")
                    elif isinstance(v, float):
                        row_vals.append(f"{v:.6g}")
                    else:
                        row_vals.append(str(v)[:40])
                # Only print row if it has any non-empty cell
                if any(rv != "" for rv in row_vals):
                    print(f"  R{row_idx}: " + " | ".join(row_vals))
    except Exception as e:
        print(f"ERROR reading {fpath}: {e}")
    print()
