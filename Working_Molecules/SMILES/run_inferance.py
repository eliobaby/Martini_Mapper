from pathlib import Path
import pandas as pd
from log_p_gnn.export import LogPGNN

CKPT = Path("exported_ckpts/extended-12-13-2024/best.ckpt")  # adjust if needed
csv_file = Path("full_list.csv")                    # input CSV with a 'cgsmiles' column

model = LogPGNN(ckpt_path=CKPT)

df = pd.read_csv(csv_file)
if "cgsmiles" not in df.columns:
    raise KeyError(f"CSV must contain a 'cgsmiles' column. Found: {list(df.columns)}")

# Ensure columns exist
for col in ["OCO", "HD", "CLF", "error"]:
    if col not in df.columns:
        df[col] = pd.NA

for i, raw in df["cgsmiles"].items():
    s = str(raw).strip() if pd.notna(raw) else ""
    if not s:
        df.at[i, "error"] = "empty cgsmiles"
        continue

    # strip surrounding quotes if present
    if (s.startswith('"') and s.endswith('"')) or (s.startswith("'") and s.endswith("'")):
        s = s[1:-1].strip()

    pred = model.predict(s)
    df.at[i, "OCO"] = float(pred.get("OCO")) if pred.get("OCO") is not None else pd.NA
    df.at[i, "HD"]  = float(pred.get("HD"))  if pred.get("HD")  is not None else pd.NA
    df.at[i, "CLF"] = float(pred.get("CLF")) if pred.get("CLF") is not None else pd.NA
    df.at[i, "error"] = ""

out_path = csv_file.with_name(csv_file.stem + "_logp.csv")
df.to_csv(out_path, index=False)
print(f"Wrote {out_path}")
