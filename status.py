import pandas as pd
import re

tsv_input = "output tsv from rt.py"
gtf_file = "annotation.gtf"
tsv_output = "tsv containing status column with included/skipped"

df = pd.read_csv(tsv_input, sep="\t")

gtf_cds = []
with open(gtf_file, "r") as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        cols = line.rstrip().split('\t')
        if len(cols) < 9 or cols[2] != "CDS":
            continue
        chrom, _, _, start, end, _, strand, _, attrs = cols
        tx_match = re.search(r'transcript_id "([^"]+)"', attrs)
        if not tx_match:
            continue
        gtf_cds.append({
            "chr": chrom,
            "start": int(start),
            "end": int(end),
            "strand": strand,
            "reference_transcript": tx_match.group(1)
        })

gtf_df = pd.DataFrame(gtf_cds)

def determine_status(row):
    chr_query = row["chr"]
    tx = row["reference_transcript"]
    start = max(1, int(row["start"]))
    end = int(row["end"])
    tx_cds = (
        gtf_df
        .loc[
            (gtf_df["chr"] == chr_query) &
            (gtf_df["reference_transcript"] == tx),
            ["start","end"]
        ]
        .sort_values("start")
    )
    if tx_cds.empty:
        return "skipped"
    for _, cds in tx_cds.iterrows():
        if start >= cds.start and end <= cds.end:
            return "included"
    return "skipped"

df["status"] = df.apply(determine_status, axis=1)
df.to_csv(tsv_output, sep="\t", index=False)
print(f"✅ Output saved to {tsv_output}")
