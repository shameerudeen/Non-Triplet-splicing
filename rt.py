#!/usr/bin/env python3
import csv
from collections import defaultdict

input_file = "tsv file containing Gene name, start and end coordinates of splicing events from JUM analysis”
annotation_file = "annotation.gtf"
output_file = "tsv with reference transcript column added along with all input columns"

rows = []
gene_ids = set()
with open(input_file, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    header = reader.fieldnames
    for rec in reader:
        rows.append(rec)
        gene_ids.add(rec['Gene'])

cds_bounds = defaultdict(lambda: defaultdict(lambda: [float('inf'), -float('inf')]))
with open(annotation_file) as gtf:
    for line in gtf:
        if line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if cols[2] != "CDS":
            continue
        start, end = int(cols[3]), int(cols[4])
        attrs = {}
        for attr in cols[8].split(';'):
            attr = attr.strip()
            if not attr:
                continue
            key, val = attr.split(' ', 1)
            attrs[key] = val.strip().strip('"')
        gene = attrs.get('gene_name')
        transcript = attrs.get('transcript_id')
        if gene in gene_ids and transcript:
            b = cds_bounds[gene][transcript]
            if start < b[0]:
                b[0] = start
            if end > b[1]:
                b[1] = end

reference = {}
for gene in gene_ids:
    best, span = "", -1
    for tid, (smin, emax) in cds_bounds[gene].items():
        if smin <= emax:
            length = emax - smin + 1
            if length > span:
                span, best = length, tid
    reference[gene] = best

with open(output_file, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(header + ['reference_transcript'])
    for rec in rows:
        writer.writerow([rec[col] for col in header] + [reference.get(rec['Gene'], "")])

