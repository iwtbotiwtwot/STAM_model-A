"""List all sheets in SU625.xlsx and STAM.xlsx, and show derivation-page candidates."""

import sys
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import openpyxl

files = [
    r"C:\Users\drwho\OneDrive\Documents\SU625\SU625.xlsx",
    r"C:\Users\drwho\OneDrive\Documents\SU625\STAM.xlsx",
]

for fpath in files:
    print("=" * 78)
    print(f"FILE: {fpath}")
    print("=" * 78)
    try:
        wb = openpyxl.load_workbook(fpath, data_only=False)
        print(f"Sheets: {wb.sheetnames}")
        for sheet_name in wb.sheetnames:
            ws = wb[sheet_name]
            print(f"  - {sheet_name}: {ws.max_row} rows x {ws.max_column} cols")
    except PermissionError as e:
        print(f"PERMISSION DENIED: {fpath} (probably open in Excel)")
    except Exception as e:
        print(f"ERROR: {e}")
    print()
