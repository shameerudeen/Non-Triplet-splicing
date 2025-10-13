#!/usr/bin/env python3
import os
import csv
from Bio import SeqIO, Seq
from collections import defaultdict

# --- CONFIGURATION ---
gtf_file = "annotation.gtf"
genome_fasta = "genome.fa"
mxe_input = "tsv containing coordinates of MXE from JUM"

nuc_output_dir = "simulated transcript sequence"
prot_output_dir = "protein sequence for both reference and alternate transcript"
idr_output_dir = "MXE region amino acid sequence specific to reference and alternate transcript used to predict disorder scores"


for d in [nuc_output_dir, prot_output_dir, idr_output_dir]:
    os.makedirs(d, exist_ok=True)


def extract_transcript_id(attr_field):
    """Extract transcript_id from GTF attributes."""
    for part in attr_field.split(";"):
        part = part.strip()
        if part.startswith('transcript_id'):
            return part.split('"')[1]
    return None

def simulate_mxe(cds_exons, exonA, exonB, status_A, status_B, mode="ref"):
    """Simulate mutually exclusive exon inclusion/exclusion."""
    if status_A == "skipped" and status_B == "skipped":
        include_A = (mode == "ref")
        include_B = (mode == "alt")
    else:
        include_A = (status_A == "included" if mode == "ref" else status_A == "skipped")
        include_B = (status_B == "included" if mode == "ref" else status_B == "skipped")

    new_exons = []
    for s, e in cds_exons:
        if not include_A and (s, e) == exonA:
            continue
        if not include_B and (s, e) == exonB:
            continue
        new_exons.append((s, e))

    if include_A and exonA not in new_exons:
        new_exons.append(exonA)
    if include_B and exonB not in new_exons:
        new_exons.append(exonB)

    return sorted(new_exons)

def build_cds_sequence(chrom, exons, strand, genome):
    """Build CDS sequence from exons and genome."""
    if chrom not in genome:
        raise ValueError(f"Chromosome {chrom} not found in genome FASTA.")
    parts = []
    exon_map = []
    for s, e in (sorted(exons) if strand == "+" else sorted(exons, reverse=True)):
        if s >= e or s < 0 or e > len(genome[chrom].seq):
            raise ValueError(f"Invalid exon coordinates: {s}-{e} on {chrom} (length {len(genome[chrom].seq)})")
        frag = genome[chrom].seq[s:e]
        if strand == "-":
            frag = frag.reverse_complement()
        parts.append(str(frag))
        exon_map.extend([(s, e)] * (e - s))
    return "".join(parts), exon_map

def map_amino_acids_to_exons(cds_seq, exon_map):
    """Map codons to their corresponding exons."""
    if len(cds_seq) < 3:
        return []
    codon_to_exon = []
    for i in range(0, len(cds_seq) - 2, 3):
        codon_exons = exon_map[i:i+3]
        if not codon_exons:
            continue
        exon_counts = defaultdict(int)
        for exon in codon_exons:
            exon_counts[exon] += 1
        majority_exon = max(exon_counts.items(), key=lambda x: x[1])[0]
        codon_to_exon.append(majority_exon)
    return codon_to_exon

def get_region_for_exon(codon_exon_list, target_exon):
    """Get residue range for a specific exon in the protein."""
    indices = [i for i, exon in enumerate(codon_exon_list) if exon == target_exon]
    if not indices:
        return None, None
    return indices[0] + 1, indices[-1] + 1

def translate_sequence(seq):
    """Translate nucleotide sequence to protein, checking for stop codons."""
    if len(seq) % 3 != 0:
        print(f"Warning: Sequence length {len(seq)} is not divisible by 3, may affect translation.")
    try:
        translated = Seq.Seq(seq).translate(to_stop=False)
        stop_codon_pos = translated.find('*')
        if stop_codon_pos != -1:
            print(f"Found in-frame stop codon at amino acid position {stop_codon_pos + 1}")
            return str(translated[:stop_codon_pos]), True
        return str(translated).replace("*", ""), False
    except Exception as e:
        print(f"Translation error: {e}")
        return "", False

print("Loading genome...")
try:
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    if not genome:
        raise ValueError("Genome FASTA file is empty or invalid.")
except FileNotFoundError:
    print(f"Error: Genome file {genome_fasta} not found.")
    exit(1)
except Exception as e:
    print(f"Error loading genome FASTA: {e}")
    exit(1)

print("Parsing GTF for CDS entries...")
cds_by_transcript = defaultdict(list)
tx_chrom = {}
tx_strand = {}

try:
    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            try:
                chrom = fields[0]
                start = int(fields[3]) - 1
                end = int(fields[4])
                strand = fields[6]
                tx = extract_transcript_id(fields[8])
                if not tx:
                    print(f"Warning: No transcript_id found in GTF line: {line.strip()}")
                    continue
                cds_by_transcript[tx].append((start, end))
                tx_chrom[tx] = chrom
                tx_strand[tx] = strand
            except ValueError as e:
                print(f"Warning: Invalid GTF line format: {line.strip()} ({e})")
                continue
except FileNotFoundError:
    print(f"Error: GTF file {gtf_file} not found.")
    exit(1)
except Exception as e:
    print(f"Error parsing GTF: {e}")
    exit(1)

print("Reading MXE data...")
try:
    with open(mxe_input) as f:
        reader = csv.DictReader(f, delimiter="\t")
        required_columns = {"Gene", "chr", "strand", "ref", "start_A", "end_A", "start_B", "end_B", "status_A", "status_B"}
        if not all(col in reader.fieldnames for col in required_columns):
            missing = required_columns - set(reader.fieldnames)
            print(f"Error: Missing required columns in TSV: {missing}")
            exit(1)
        mxe_rows = [row for row in reader]
    if not mxe_rows:
        print(f"Error: No data rows found in {mxe_input}.")
        exit(1)
except FileNotFoundError:
    print(f"Error: MXE file {mxe_input} not found.")
    exit(1)
except Exception as e:
    print(f"Error reading MXE TSV: {e}")
    exit(1)

gene_counts = defaultdict(int)
for row in mxe_rows:
    gene_counts[row["Gene"]] += 1

skipped_rows = 0
processed_rows = 0
for row in mxe_rows:
    try:
        gene = row["Gene"]
        chrom = row["chr"]
        strand = row["strand"]
        tx = row["ref"]
        start_A = int(row["start_A"])
        end_A = int(row["end_A"])
        start_B = int(row["start_B"])
        end_B = int(row["end_B"])
        status_A = row["status_A"].lower()
        status_B = row["status_B"].lower()

        if tx not in cds_by_transcript:
            print(f"Warning: Transcript {tx} not found in GTF for gene {gene}, skipping.")
            skipped_rows += 1
            continue
        if tx_chrom[tx] != chrom or tx_strand[tx] != strand:
            print(f"Warning: Mismatch in chrom ({tx_chrom[tx]} vs {chrom}) or strand ({tx_strand[tx]} vs {strand}) for {tx}, skipping.")
            skipped_rows += 1
            continue
        if status_A not in {"included", "skipped"} or status_B not in {"included", "skipped"}:
            print(f"Warning: Invalid status_A ({status_A}) or status_B ({status_B}) for {gene}, skipping.")
            skipped_rows += 1
            continue
        if start_A >= end_A or start_B >= end_B:
            print(f"Warning: Invalid exon coordinates for {gene}: A({start_A}-{end_A}), B({start_B}-{end_B}), skipping.")
            skipped_rows += 1
            continue

        exonA = (start_A - 1, end_A)
        exonB = (start_B - 1, end_B)
        len_A = end_A - start_A + 1
        len_B = end_B - start_B + 1

        exclude_flag = False
        region_data = {}
        cds_seqs = {}
        protein_seqs = {}

        for mode, exon_included in [("ref", exonA if status_A == "included" else exonB),
                                    ("alt", exonB if status_A == "included" else exonA)]:
            try:
                cds_exons = simulate_mxe(cds_by_transcript[tx], exonA, exonB, status_A, status_B, mode)
                if not cds_exons:
                    print(f"Warning: No exons generated for {gene}|{mode}, skipping row.")
                    exclude_flag = True
                    break
                cds_seq, exon_map = build_cds_sequence(chrom, cds_exons, strand, genome)
                if not cds_seq:
                    print(f"Warning: Empty CDS sequence for {gene}|{mode}, skipping row.")
                    exclude_flag = True
                    break

                prot, has_stop = translate_sequence(cds_seq)
                if not prot:
                    print(f"Warning: Empty protein sequence for {gene}|{mode}, skipping row.")
                    exclude_flag = True
                    break
                codon_exon_map = map_amino_acids_to_exons(cds_seq, exon_map)
                res_start, res_end = get_region_for_exon(codon_exon_map, exon_included)
                region_seq = prot[res_start - 1:res_end] if res_start else ""

                exon_nt_len = len_A if exon_included == exonA else len_B
                aa_codon_len = len(region_seq) * 3
                diff = abs(aa_codon_len - exon_nt_len)
                is_triplet = exon_nt_len % 3 == 0
                exon_seq = genome[chrom].seq[exon_included[0]:exon_included[1]]
                if strand == "-":
                    exon_seq = exon_seq.reverse_complement()
                exon_seq = str(exon_seq)
                if exon_nt_len >= 150 and diff > 6:
                    print(f"[FLAG] {gene}|{mode}: exon {exon_included} is {exon_nt_len} nt, but translated region is only {len(region_seq)} aa")
                    print(f"  - Exon sequence: {exon_seq[:50]}... (first 50 nt)")
                    print(f"  - Is triplet: {is_triplet} (length {exon_nt_len} nt)")
                    print(f"  - In-frame stop codon: {'Yes' if has_stop else 'No'}")
                    exclude_flag = True

                region_data[mode] = (res_start, res_end, region_seq)
                cds_seqs[mode] = cds_seq
                protein_seqs[mode] = prot
            except Exception as e:
                print(f"Error processing {gene}|{mode}: {e}")
                exclude_flag = True
                break

        if exclude_flag:
            print(f"Skipping FASTA output for {gene} due to length mismatch or processing error.")
            skipped_rows += 1
            continue

        safe_gene = "".join(c if c.isalnum() or c in "-_" else "_" for c in gene)
        fname = f"{safe_gene}_{start_A}_{end_A}_{start_B}_{end_B}.fa" if gene_counts[gene] > 1 else f"{safe_gene}.fa"

        with open(os.path.join(nuc_output_dir, fname), "w") as f:
            for mode in ["ref", "alt"]:
                f.write(f">{safe_gene}|{mode}\n")
                seq = cds_seqs[mode]
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")

        with open(os.path.join(prot_output_dir, fname), "w") as f:
            for mode in ["ref", "alt"]:
                f.write(f">{safe_gene}|{mode}\n")
                seq = protein_seqs[mode]
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")

        with open(os.path.join(idr_output_dir, fname), "w") as f:
            for mode in ["ref", "alt"]:
                res_start, res_end, region_seq = region_data[mode]
                f.write(f">{safe_gene}|{mode}\n")
                f.write(f"Mutually exclusive region: residues {res_start} to {res_end} (length {res_end - res_start + 1 if res_start else 0})\n")
                for i in range(0, len(region_seq), 60):
                    f.write(region_seq[i:i+60] + "\n")

        processed_rows += 1

    except Exception as e:
        print(f"Error processing row for gene {row.get('Gene', 'unknown')}: {e}")
        skipped_rows += 1
        continue

print(f"\n✅ All done. Processed {processed_rows} rows, skipped {skipped_rows} rows due to flags or errors.")

