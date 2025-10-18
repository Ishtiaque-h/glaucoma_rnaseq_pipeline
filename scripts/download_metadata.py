# scripts/download_metadata.py
import os
import pandas as pd
from pysradb.sraweb import SRAweb

def fetch_metadata(SRA_ID: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    db = SRAweb()
    df = db.sra_metadata(SRA_ID, assay_type="RNA-Seq", detailed=True)
    output_path = os.path.join(out_dir, f"{SRA_ID}_metadata.csv")
    df.to_csv(output_path, index=False)
    print(f"Metadata saved to {output_path}")
    return df

if __name__ == "__main__":
    SRA_ID = "SRP394552"
    OUTPUT_DIR = "data/metadata"
    fetch_metadata(SRA_ID, OUTPUT_DIR)
