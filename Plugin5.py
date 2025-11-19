import os
import pandas as pd

# Source and destination directories
src_dir  = r"C:\Users\WerlandL\FRET2"
dest_dir = r"C:\Users\WerlandL\FRET3"
os.makedirs(dest_dir, exist_ok=True)

# File name components
prefixes = ["Me1xy", "Me2xy", "Me3xy", "Me4xy", "Me5xy", "Me6xy", "Me7xy", "Me8xy"]
types    = ["ResultsIN-Ble", "ResultsOUT"]
numbers  = range(1, 21)  # 01 … 20

# Helper to convert pandas‐col‐index → Excel letter
def col_letter(idx):
    idx += 1
    s = ""
    while idx:
        idx, rem = divmod(idx - 1, 26)
        s = chr(65 + rem) + s
    return s

# Helper to slice Excel rows in pandas
def excel_slice(df, col, start, end):
    # header row = 1 → iloc index 0, so Excel N → iloc N-2
    return df[col].iloc[start-2 : end-1]

# The four averaging ranges and the labels for A59–A62
avg_ranges = [(26,28), (29,31), (2,3), (54,55)]
labels1    = ["DonBB", "DonAB", "AccBB", "AccAB"]

# The labels for A63–A72 and where their formulas/constants go
labels2 = [
    "BG don",    # row 63 → 170
    "BG acc",    # row 64 → 135
    "DonBB-BG",  # row 65 → =X59-X63
    "DonAB-BG",  # row 66 → =X60-X63
    "AccBB-BG",  # row 67 → =X61-X64
    "ACCAB-BG",  # row 68 → =X62-X64
    "FRET",      # row 69 → =100*(X66-X65)/X66
    "Bleached",  # row 70 → =100*(X67-X68)/X67
    "Fac",       # row 71 → =(100-X70)/100
    "CorrFRET"   # row 72 → =100*(X66-X65)/(X66-X71*X65)
]

for prefix in prefixes:
    for num in numbers:
        for typ in types:
            orig_fname = f"{prefix}{num:02d}{typ}.csv"
            orig_path  = os.path.join(src_dir, orig_fname)
            if not os.path.exists(orig_path):
                print(f"→ Skipping missing: {orig_path}")
                continue

            # 1) Load data
            df = pd.read_csv(orig_path)
            target_cols = df.columns[1::2]

            # 2) Compute the four average‐rows (→ Excel 59–62)
            avg_dicts = []
            for start, end in avg_ranges:
                d = {c: excel_slice(df, c, start, end).mean() for c in target_cols}
                avg_dicts.append(d)
            new_rows_df = pd.DataFrame(avg_dicts, columns=df.columns)
            new_rows_df[df.columns[0]] = labels1  # A59–A62

            # 3) Determine Excel row numbers for the appended block
            n_orig = len(df)
            first  = n_orig + 2  # new_rows_df.iloc[0] → Excel row 59
            r59, r60, r61, r62 = first, first+1, first+2, first+3
            r63, r64 = first+4, first+5
            r65, r66, r67, r68 = first+6, first+7, first+8, first+9
            r69, r70, r71, r72 = first+10, first+11, first+12, first+13

            # 4) Build rows 63–72 with labels2 in col A and formulas/constants
            label_rows = []
            for lbl in labels2:
                row = {c: pd.NA for c in df.columns}
                row[df.columns[0]] = lbl
                for c in target_cols:
                    L = col_letter(df.columns.get_loc(c))
                    if   lbl == "BG don":
                        row[c] = "170"
                    elif lbl == "BG acc":
                        row[c] = "135"
                    elif lbl == "DonBB-BG":
                        row[c] = f"={L}{r59}-{L}{r63}"
                    elif lbl == "DonAB-BG":
                        row[c] = f"={L}{r60}-{L}{r63}"
                    elif lbl == "AccBB-BG":
                        row[c] = f"={L}{r61}-{L}{r64}"
                    elif lbl == "ACCAB-BG":
                        row[c] = f"={L}{r62}-{L}{r64}"
                    elif lbl == "FRET":
                        row[c] = f"=100*({L}{r66}-{L}{r65})/{L}{r66}"
                    elif lbl == "Bleached":
                        row[c] = f"=100*({L}{r67}-{L}{r68})/{L}{r67}"
                    elif lbl == "Fac":
                        row[c] = f"=(100-{L}{r70})/100"
                    elif lbl == "CorrFRET":
                        row[c] = f"=100*({L}{r66}-{L}{r65})/({L}{r66}-{L}{r71}*{L}{r65})"
                label_rows.append(row)
            label_rows_df = pd.DataFrame(label_rows, columns=df.columns)

            # 5) Concatenate and save
            df_out = pd.concat([df, new_rows_df, label_rows_df], ignore_index=True)
            new_fname = "Ed_" + orig_fname
            new_path  = os.path.join(dest_dir, new_fname)
            df_out.to_csv(new_path, index=False)
            print(f"✅ Processed: {new_path}")
