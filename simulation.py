import os
import csv
from collections import Counter
from copy import deepcopy
from Bio import SeqIO

# CONFIG
tsv_input = "output tsv from status.py"
gtf_file = "annotation.gtf"
genome_fasta = "genome.fa"
tsv_output = "tsv with all parameters described in the methods"
outcds_tsv = "tsv for splicing events outside of the cds bpundary"
false_exons_tsv = "tsv for more than one annotated exon within JUM start and end coordinates"
error_tsv = "tsv for excluded splicing events"
sim_seq_dir = "directory containing simulated alternative transcript sequence along with the region simulated"


os.makedirs(sim_seq_dir, exist_ok=True)


def extract_transcript_id(attr_field):
    for part in attr_field.split(";"):
        part = part.strip()
        if part.startswith('transcript_id'):
            return part.split('"')[1]
    return None

def is_within_cds_bounds(rs, re, cds_exons):
    if not cds_exons:
        return False, "no CDS on transcript"
    min_s = min(s for s, e in cds_exons)
    max_e = max(e for s, e in cds_exons)
    if rs < min_s or re > max_e:
        return False, "region outside CDS bounds"
    return True, "within CDS bounds"

def find_false_exons(rs, re, cds_exons):
    overlapping_exons = [(s, e) for s, e in cds_exons if s <= re and e > rs]
    return len(overlapping_exons) > 1, overlapping_exons

def simulate_region(exons, rs, re, action):
    new_exons = []
    if action == "exclude":
        for s, e in exons:
            if re <= s or rs >= e:
                new_exons.append((s, e))
            else:
                if rs > s:
                    new_exons.append((s, min(rs, e)))
                if re < e:
                    new_exons.append((max(re, s), e))
    else:
        new_exons = exons[:] + [(rs, re)]
    return new_exons

def find_stop(full_seq, cds_exons, strand, utr_tail_length):
    if not cds_exons:
        return "No CDS exons"
    seq = full_seq.upper()
    pos_atg = seq.find("ATG")
    if pos_atg < 0:
        return "No start codon"
    stops = {"TAA", "TAG", "TGA"}
    stop_i = None
    for i in range(pos_atg, len(seq) - 2, 3):
        if seq[i:i+3] in stops:
            stop_i = i
            break
    if stop_i is None:
        return "No stop codon"
    cur = 0
    for s, e in cds_exons:
        length = e - s
        if cur <= stop_i < cur + length:
            off = stop_i - cur
            if strand == "+":
                return str(s + off + 1)  # Convert to 1-based
            else:
                return str(e - off)  # Convert to 1-based
        cur += length
    if cur <= stop_i < cur + utr_tail_length:
        off = stop_i - cur
        if strand == "+":
            last_cds_end = cds_exons[-1][1]
            return str(last_cds_end + off + 1)  # Convert to 1-based
        else:
            last_cds_start = cds_exons[-1][0]
            return str(last_cds_start - off)  # Convert to 1-based
    return "Mapping error"

def calculate_utr_distance(sc_pos, exons, strand, end_utr, full_exons, chrom, genome):
    try:
        sc_pos = int(sc_pos)
    except ValueError:
        return 0
    sorted_exons = sorted(exons, reverse=(strand == "-"))
    sorted_full_exons = sorted(full_exons, reverse=(strand == "-"))
    parts = []
    utr_tail_length = 0
    for s, e in sorted_exons:
        frag = genome[chrom].seq[s:e]
        if strand == "-":
            frag = frag.reverse_complement()
        parts.append(str(frag))
    if strand == "+":
        last_cds_end = sorted_exons[-1][1]
        transcript_end = end_utr
        if transcript_end > last_cds_end:
            tail = genome[chrom].seq[last_cds_end:transcript_end]
            parts.append(str(tail))
            utr_tail_length = len(tail)
    else:
        last_cds_start = sorted_exons[-1][0]
        transcript_end = end_utr
        if last_cds_start > transcript_end:
            tail = genome[chrom].seq[transcript_end:last_cds_start].reverse_complement()
            parts.append(str(tail))
            utr_tail_length = len(tail)
    full_seq = "".join(parts)
    cur = 0
    stop_i = None
    for s, e in sorted_exons:
        length = e - s
        if strand == "+":
            if s <= sc_pos - 1 < e:  # Adjust for 1-based input
                stop_i = cur + (sc_pos - 1 - s)
                break
        else:
            if s <= sc_pos - 1 < e:  # Adjust for 1-based input
                stop_i = cur + (e - (sc_pos - 1) - 1)
                break
        cur += length
    if stop_i is None:
        if strand == "+":
            if last_cds_end <= sc_pos - 1 <= transcript_end:
                stop_i = cur + (sc_pos - 1 - last_cds_end)
        else:
            if transcript_end <= sc_pos - 1 <= last_cds_start:
                stop_i = cur + (last_cds_start - (sc_pos - 1))
    if stop_i is None:
        return 0
    distance = len(full_seq) - stop_i
    if strand == "+":
        if sc_pos - 1 > end_utr:
            distance = -distance
    else:
        if sc_pos - 1 < end_utr:
            distance = -distance
    return distance

def calculate_transcript_distance(sc_pos, exons, strand, junc, full_exons, chrom, genome):
    try:
        sc_pos = int(sc_pos)
    except ValueError:
        return 0
    sorted_exons = sorted(exons, reverse=(strand == "-"))
    sorted_full_exons = sorted(full_exons, reverse=(strand == "-"))
    parts = []
    utr_tail_length = 0
    for s, e in sorted_exons:
        frag = genome[chrom].seq[s:e]
        if strand == "-":
            frag = frag.reverse_complement()
        parts.append(str(frag))
    if strand == "+":
        last_cds_end = sorted_exons[-1][1]
        transcript_end = max(e for s, e in sorted_full_exons)
        if transcript_end > last_cds_end:
            tail = genome[chrom].seq[last_cds_end:transcript_end]
            parts.append(str(tail))
            utr_tail_length = len(tail)
    else:
        last_cds_start = sorted_exons[-1][0]
        transcript_end = min(s for s, e in sorted_full_exons)
        if last_cds_start > transcript_end:
            tail = genome[chrom].seq[transcript_end:last_cds_start].reverse_complement()
            parts.append(str(tail))
            utr_tail_length = len(tail)
    full_seq = "".join(parts)
    cur = 0
    sc_i = None
    for s, e in sorted_exons:
        length = e - s
        if strand == "+":
            if s <= sc_pos - 1 < e:  # Adjust for 1-based input
                sc_i = cur + (sc_pos - 1 - s)
                break
        else:
            if s <= sc_pos - 1 < e:  # Adjust for 1-based input
                sc_i = cur + (e - (sc_pos - 1) - 1)
                break
        cur += length
    if sc_i is None:
        if strand == "+":
            if last_cds_end <= sc_pos - 1 <= transcript_end:
                sc_i = cur + (sc_pos - 1 - last_cds_end)
        else:
            if transcript_end <= sc_pos - 1 <= last_cds_start:
                sc_i = cur + (last_cds_start - (sc_pos - 1))
    if sc_i is None:
        return 0
    junc_i = None
    cur = 0
    for s, e in sorted_exons:
        length = e - s
        if strand == "+":
            if s <= junc - 1 < e:  # Adjust for 1-based input
                junc_i = cur + (junc - 1 - s)
                break
        else:
            if s <= junc - 1 < e:  # Adjust for 1-based input
                junc_i = cur + (e - (junc - 1) - 1)
                break
        cur += length
    if junc_i is None:
        if strand == "+":
            if last_cds_end <= junc - 1 <= transcript_end:
                junc_i = cur + (junc - 1 - last_cds_end)
        else:
            if transcript_end <= junc - 1 <= last_cds_start:
                junc_i = cur + (last_cds_start - (junc - 1))
    if junc_i is None:
        return 0
    distance = junc_i - sc_i
    
    if (strand == "+" and sc_pos - 1 < junc - 1) or (strand == "-" and sc_pos - 1 > junc - 1):
        distance = abs(distance)  # Upstream: positive
    else:
        distance = -abs(distance)  # Downstream: negative
    return distance


print("Loading genome…")
try:
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
except FileNotFoundError:
    raise FileNotFoundError(f"Genome FASTA file not found: {genome_fasta}")
except Exception as e:
    raise Exception(f"Error parsing genome FASTA: {str(e)}")


print("Parsing GTF (exons + CDS)…")
gtf_exons = {}
transcript_cds = {}
transcript_to_id = {}
try:
    with open(gtf_file) as gtf:
        for L in gtf:
            if L.startswith("#"):
                continue
            f = L.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, typ = f[0], f[2]
            s0, e0 = int(f[3]) - 1, int(f[4])
            strand_gtf = f[6]
            txid = extract_transcript_id(f[8])
            if not txid:
                continue
            if typ == "exon":
                gtf_exons.setdefault(txid, {"chrom": chrom, "strand": strand_gtf, "exons": []})
                gtf_exons[txid]["exons"].append((s0, e0))
            elif typ == "CDS":
                transcript_cds.setdefault(txid, []).append((s0, e0))
except FileNotFoundError:
    raise FileNotFoundError(f"GTF file not found: {gtf_file}")
except Exception as e:
    raise Exception(f"Error parsing GTF file: {str(e)}")

print(f"  → {len(gtf_exons)} transcripts with exons")
print(f"  → {len(transcript_cds)} transcripts with CDS")


print("Detecting gene column and validating header…")
try:
    with open(tsv_input) as inf:
        header = inf.readline().rstrip("\n").split("\t")
except FileNotFoundError:
    raise FileNotFoundError(f"Input TSV file not found: {tsv_input}")
except Exception as e:
    raise Exception(f"Error reading input TSV: {str(e)}")


required_columns = ["start", "end", "strand", "reference_transcript", "status", "Group"]
optional_columns = ["chr", "ID", "avgwtPSI", "avgsmgPSI", "deltaPSI"]
possible_gene_cols = ["Gene", "gene", "gene_name", "gene_id"]
gene_col = None
for col in possible_gene_cols:
    if col in header:
        gene_col = col
        break
else:
    raise KeyError(f"No gene column found; header is: {header!r}")

missing_required = [col for col in required_columns if col not in header]
if missing_required:
    raise KeyError(f"Missing required columns: {missing_required}; header is: {header!r}")


main_fields = [
    gene_col, "ID", "chr", "strand", "start", "end", "In_frame_sc", "avgwtPSI", "avgsmgPSI",
    "deltaPSI", "Group", "reference_transcript", "status", "sc_ref", "sc_alt",
    "exon_after_sc_alt", "exon_after_end", "recoded_exons", "ntd_recoded_ref",
    "ntds_altv_rec", "dual_coding", "exon_junc", "distance_sc_alt", "distance_sc_ref",
    "end_UTR", "dist_ref_UTR", "dist_alt_UTR", "length", "category", "Abs_deltaPSI",
    "Phase", "AS_event_ID", "NMD", "Significance"
]

false_exons_fields = header + ["exclusion_reason", "overlapping_exons"]
error_fields = main_fields + ["error_reason"]

try:
    with open(tsv_input) as inf:
        reader = csv.DictReader(inf, delimiter="\t")
        rdr_rows = list(reader)
except Exception as e:
    raise Exception(f"Error parsing input TSV: {str(e)}")

if not rdr_rows:
    print("Warning: Input TSV is empty (no data rows). Creating empty output files.")

genes = [r[gene_col] for r in rdr_rows]
counts = Counter(genes)


for row in rdr_rows:
    txid = row["reference_transcript"]
    id_val = row.get("ID", "")
    if id_val and txid:
        transcript_to_id[txid] = id_val


print("Processing rows…")
false_exons_results = []
error_results = []
try:
    with open(tsv_output, "w", newline="") as outf, \
         open(outcds_tsv, "w", newline="") as outcds, \
         open(error_tsv, "w", newline="") as errorf:
        base_fields = header
        cds_fields = base_fields + ["exclusion_reason"]
        w_main = csv.DictWriter(outf, fieldnames=main_fields, delimiter="\t")
        w_excluded = csv.DictWriter(outcds, fieldnames=cds_fields, delimiter="\t")
        w_error = csv.DictWriter(errorf, fieldnames=error_fields, delimiter="\t")
        w_main.writeheader()
        w_excluded.writeheader()
        w_error.writeheader()

        for row in rdr_rows:
            error_reason = None
            tx = row["reference_transcript"]
            gene = row[gene_col]
            try:
                start = int(row["start"])
                end = int(row["end"])
            except ValueError:
                error_reason = "Invalid start/end coordinates"
                row.update({f: row.get(f, "") for f in main_fields if f not in row})
                row.update({"exon_after_end": "", "sc_ref": "invalid coordinates", "sc_alt": "",
                            "recoded_exons": "", "ntd_recoded_ref": "", "ntds_altv_rec": "",
                            "dual_coding": "", "exon_after_sc_alt": "", "exon_junc": "",
                            "distance_sc_alt": "", "distance_sc_ref": "", "In_frame_sc": "",
                            "end_UTR": "", "dist_ref_UTR": "", "dist_alt_UTR": ""})
                row["error_reason"] = error_reason
                error_results.append(row)
                continue
            rs, re = start - 1, end
            status = row["status"].lower()
            strand = row["strand"]

           
            if strand not in ["+", "-"]:
                error_reason = f"Invalid strand value '{strand}'"
                row.update({f: row.get(f, "") for f in main_fields if f not in row})
                row.update({"exon_after_end": "", "sc_ref": "invalid strand", "sc_alt": "",
                            "recoded_exons": "", "ntd_recoded_ref": "", "ntds_altv_rec": "",
                            "dual_coding": "", "exon_after_sc_alt": "", "exon_junc": "",
                            "distance_sc_alt": "", "distance_sc_ref": "", "In_frame_sc": "",
                            "end_UTR": "", "dist_ref_UTR": "", "dist_alt_UTR": ""})
                row["error_reason"] = error_reason
                error_results.append(row)
                continue

            info = gtf_exons.get(tx)
            if not info:
                error_reason = "No exons for transcript"
                row.update({f: row.get(f, "") for f in main_fields if f not in row})
                row.update({"exon_after_end": "", "sc_ref": "no exons", "sc_alt": "",
                            "recoded_exons": "", "ntd_recoded_ref": "", "ntds_altv_rec": "",
                            "dual_coding": "", "exon_after_sc_alt": "", "exon_junc": "",
                            "distance_sc_alt": "", "distance_sc_ref": "", "In_frame_sc": "",
                            "end_UTR": "", "dist_ref_UTR": "", "dist_alt_UTR": ""})
                row["error_reason"] = error_reason
                error_results.append(row)
                continue

            chrom = info["chrom"]
            cds_exons = transcript_cds.get(tx, [])
            full_exons = info["exons"]
            if not cds_exons:
                error_reason = "No CDS for transcript"
                row.update({f: row.get(f, "") for f in main_fields if f not in row})
                row.update({"exon_after_end": "", "sc_ref": "no CDS", "sc_alt": "",
                            "recoded_exons": "", "ntd_recoded_ref": "", "ntds_altv_rec": "",
                            "dual_coding": "", "exon_after_sc_alt": "", "exon_junc": "",
                            "distance_sc_alt": "", "distance_sc_ref": "", "In_frame_sc": "",
                            "end_UTR": "", "dist_ref_UTR": "", "dist_alt_UTR": ""})
                row["error_reason"] = error_reason
                error_results.append(row)
                continue

            
            is_false, overlapping_exons = find_false_exons(rs, re, cds_exons)
            if is_false:
                excluded = {c: row[c] for c in base_fields}
                excluded["exclusion_reason"] = "multiple CDS exons overlap"
                excluded["overlapping_exons"] = ";".join([f"{s}-{e}" for s, e in overlapping_exons])
                false_exons_results.append(excluded)
                continue

            
            within, reason = is_within_cds_bounds(rs, re, cds_exons)
            if not within:
                excluded = {c: row[c] for c in base_fields}
                excluded["exclusion_reason"] = reason
                w_excluded.writerow(excluded)
                continue

            
            exon_after_end = 0
            ref_pos = start if strand == "-" else end
            if strand == "+":
                for s, e in full_exons:
                    if s > ref_pos:
                        exon_after_end += 1
            else:
                for s, e in full_exons:
                    if e - 1 < ref_pos:
                        exon_after_end += 1
            row["exon_after_end"] = str(exon_after_end)

            # prepare sorted lists
            cds_sorted = sorted(cds_exons, reverse=(strand == "-"))
            full_sorted = sorted(full_exons, reverse=(strand == "-"))

           
            parts = []
            utr_tail_ref = 0
            for s, e in cds_sorted:
                frag = genome[chrom].seq[s:e]
                if strand == "-":
                    frag = frag.reverse_complement()
                parts.append(str(frag))
            if strand == "+":
                last_end = cds_sorted[-1][1]
                if full_sorted[-1][1] > last_end:
                    tail = genome[chrom].seq[last_end:full_sorted[-1][1]]
                    parts.append(str(tail))
                    utr_tail_ref = len(tail)
            else:
                first_start = full_sorted[-1][0]
                last_start = cds_sorted[-1][0]
                if last_start > first_start:
                    tail = genome[chrom].seq[first_start:last_start].reverse_complement()
                    parts.append(str(tail))
                    utr_tail_ref = len(tail)
            ref_seq = "".join(parts)
            row["sc_ref"] = find_stop(ref_seq, cds_sorted, strand, utr_tail_ref)

           
            action = "include" if status == "skipped" else "exclude"
            alt_exons = simulate_region(deepcopy(cds_exons), rs, re, action)
            alt_sorted = sorted(alt_exons, reverse=(strand == "-"))
            parts = []
            utr_tail_alt = 0
            for s, e in alt_sorted:
                frag = genome[chrom].seq[s:e]
                if strand == "-":
                    frag = frag.reverse_complement()
                parts.append(str(frag))
            if strand == "+":
                last_end = alt_sorted[-1][1]
                if full_sorted[-1][1] > last_end:
                    tail = genome[chrom].seq[last_end:full_sorted[-1][1]]
                    parts.append(str(tail))
                    utr_tail_alt = len(tail)
            else:
                first_start = full_sorted[-1][0]
                last_start = alt_sorted[-1][0]
                if last_start > first_start:
                    tail = genome[chrom].seq[first_start:last_start].reverse_complement()
                    parts.append(str(tail))
                    utr_tail_alt = len(tail)
            alt_seq = "".join(parts)
            row["sc_alt"] = find_stop(alt_seq, alt_sorted, strand, utr_tail_alt)

            # Write FASTA for simulated sequence
            if alt_seq:
                safe_gene = "".join(c if c.isalnum() or c in "-_" else "_" for c in gene)
                if counts[gene] > 1:
                    fasta_filename = f"{safe_gene}_{start}_{end}.fa"
                else:
                    fasta_filename = f"{safe_gene}.fa"
                fasta_path = os.path.join(sim_seq_dir, fasta_filename)
                try:
                    with open(fasta_path, "w") as fasta_out:
                        # Write simulated region sequence
                        sim_region_seq = genome[chrom].seq[rs:re]
                        if strand == "-":
                            sim_region_seq = sim_region_seq.reverse_complement()
                        sim_region_seq = str(sim_region_seq)
                        fasta_header_region = f">{safe_gene}|tx:{tx}|chr:{chrom}|strand:{strand}|action:{action}|simulated_region"
                        fasta_out.write(f"{fasta_header_region}\n")
                        for i in range(0, len(sim_region_seq), 60):
                            fasta_out.write(f"{sim_region_seq[i:i+60]}\n")
                        # Write simulated region length
                        fasta_header_length = f">{safe_gene}|tx:{tx}|chr:{chrom}|strand:{strand}|action:{action}|region_length"
                        fasta_out.write(f"{fasta_header_length}\n")
                        fasta_out.write(f"{len(sim_region_seq)}\n")
                        # Write complete alternate sequence
                        fasta_header_alt = f">{safe_gene}|tx:{tx}|chr:{chrom}|strand:{strand}|action:{action}|alternate_sequence"
                        fasta_out.write(f"{fasta_header_alt}\n")
                        for i in range(0, len(alt_seq), 60):
                            fasta_out.write(f"{alt_seq[i:i+60]}\n")
                        # Write complete reference sequence
                        fasta_header_ref = f">{safe_gene}|tx:{tx}|chr:{chrom}|strand:{strand}|action:{action}|reference_sequence"
                        fasta_out.write(f"{fasta_header_ref}\n")
                        for i in range(0, len(ref_seq), 60):
                            fasta_out.write(f"{ref_seq[i:i+60]}\n")
                except Exception as e:
                    print(f"Warning: Failed to write FASTA {fasta_path}: {str(e)}")

           
            try:
                sc_alt = int(row["sc_alt"])
                row["In_frame_sc"] = "YES" if start <= sc_alt <= end else "NO"
            except ValueError:
                row["In_frame_sc"] = "NO"

           
            recoded_exons = 0
            try:
                sc_pos = int(row["sc_alt"])
                ref_pos = start if strand == "-" else end
                if strand == "+":
                    for s, e in cds_exons:
                        if ref_pos < s <= sc_pos:
                            recoded_exons += 1
                        elif s <= sc_pos <= e:
                            recoded_exons += 1
                else:
                    for s, e in cds_exons:
                        if e - 1 < ref_pos and e - 1 >= sc_pos:
                            recoded_exons += 1
                        elif s <= sc_pos <= e:
                            recoded_exons += 1
            except ValueError:
                recoded_exons = 0
            row["recoded_exons"] = str(recoded_exons)

            
            ntd_recoded_ref = 0
            try:
                sc_pos = int(row["sc_ref"])
                exons_to_use = cds_exons if status == "included" else alt_exons
                if strand == "+":
                    if sc_pos > end:
                        for s, e in exons_to_use:
                            if end < s < sc_pos and not (s <= sc_pos <= e):
                                ntd_recoded_ref += e - s
                            elif s <= sc_pos <= e:
                                ntd_recoded_ref += sc_pos - max(s, end)
                    else:
                        for s, e in exons_to_use:
                            if sc_pos < s < end and not (s <= sc_pos <= e):
                                ntd_recoded_ref += e - s
                            elif s <= sc_pos <= e:
                                ntd_recoded_ref += min(e, end) - sc_pos
                else:
                    if sc_pos < start:
                        for s, e in exons_to_use:
                            if sc_pos < e <= start and not (s <= sc_pos <= e):
                                ntd_recoded_ref += e - s
                            elif s <= sc_pos <= e:
                                ntd_recoded_ref += min(e, start) - sc_pos
                    else:
                        for s, e in exons_to_use:
                            if start < e <= sc_pos and not (s <= sc_pos <= e):
                                ntd_recoded_ref += e - s
                            elif s <= sc_pos <= e:
                                ntd_recoded_ref += sc_pos - max(s, start)
            except ValueError:
                ntd_recoded_ref = 0
            row["ntd_recoded_ref"] = str(ntd_recoded_ref)

           
            ntds_altv_rec = 0
            try:
                sc_pos = int(row["sc_alt"])
                exons_to_use = cds_exons if status == "included" else alt_exons
                if strand == "+":
                    if sc_pos > start:
                        for s, e in exons_to_use:
                            if s <= sc_pos and e > start:
                                segment_start = max(s, start)
                                segment_end = min(e, sc_pos)
                                if segment_end > segment_start:
                                    ntds_altv_rec += segment_end - segment_start
                    else:
                        for s, e in exons_to_use:
                            if sc_pos < s < start and not (s <= sc_pos <= e):
                                ntds_altv_rec += e - s
                            elif s <= sc_pos <= e:
                                ntds_altv_rec += min(e, start) - sc_pos
                else:
                    if sc_pos < end:
                        for s, e in exons_to_use:
                            if s <= end and e > sc_pos:
                                segment_start = max(s, sc_pos)
                                segment_end = min(e, end)
                                if segment_end > segment_start:
                                    ntds_altv_rec += segment_end - segment_start
                    else:
                        for s, e in exons_to_use:
                            if s <= sc_pos and e > end:
                                segment_start = max(s, end)
                                segment_end = min(e, sc_pos)
                                if segment_end > segment_start:
                                    ntds_altv_rec += segment_end - segment_start
            except ValueError:
                ntds_altv_rec = 0
            row["ntds_altv_rec"] = str(ntds_altv_rec)

            
            dual_coding = 0
            try:
                sc_pos = int(row["sc_alt"])
                if start <= sc_pos <= end or row["Group"] == "Triplets":
                    dual_coding = 0
                else:
                    exons_to_use = cds_exons if status == "included" else alt_exons
                    if strand == "+":
                        if sc_pos > end:
                            for s, e in exons_to_use:
                                if end < s < sc_pos and not (s <= sc_pos <= e):
                                    dual_coding += e - s
                                elif s <= sc_pos <= e:
                                    dual_coding += sc_pos - max(s, end)
                        else:
                            for s, e in exons_to_use:
                                if sc_pos < s < end and not (s <= sc_pos <= e):
                                    dual_coding += e - s
                                elif s <= sc_pos <= e:
                                    dual_coding += min(e, end) - sc_pos
                    else:
                        if sc_pos < start:
                            for s, e in exons_to_use:
                                if sc_pos < e <= start and not (s <= sc_pos <= e):
                                    dual_coding += e - s
                                elif s <= sc_pos <= e:
                                    dual_coding += min(e, start) - sc_pos
                        else:
                            for s, e in exons_to_use:
                                if start < e <= sc_pos and not (s <= sc_pos <= e):
                                    dual_coding += e - s
                                elif s <= sc_pos <= e:
                                    dual_coding += sc_pos - max(s, start)
            except ValueError:
                dual_coding = 0
            row["dual_coding"] = str(dual_coding)

            
            downstream = 0
            try:
                sc_pos = int(row["sc_alt"])
                if strand == "+":
                    for s, e in full_exons:
                        if s > sc_pos:
                            downstream += 1
                else:
                    for s, e in full_exons:
                        if e - 1 < sc_pos:
                            downstream += 1
            except ValueError:
                pass
            row["exon_after_sc_alt"] = str(downstream)

           
            if row.get("category") == "a3ss":
                if len(full_exons) >= 2:
                    if strand == "+":
                        # Sort by exon end position; take the second-to-last exon
                        penultimate_end = sorted(full_exons, key=lambda x: x[1])[-2][1]
                        junc = penultimate_end
                    else:
                        # Sort by exon start (reverse order); take second-to-last exon start
                        penultimate_start = sorted(full_exons, key=lambda x: x[0], reverse=True)[-2][0]
                        junc = penultimate_start + 1  # FIX: add 1 to represent the correct junction
                else:
                    if strand == "+":
                        # If only one exon, junction is right after the last exon end
                        junc = min(s for s, e in full_exons if e == max(e2 for s2, e2 in full_exons)) + 1
                    else:
                        # FIX: Add 1 to correct junction position for - strand
                        junc = max(e for s, e in full_exons if s == min(s2 for s2, e2 in full_exons)) + 1
            else:
                if strand == "+":
                    junc = min(s for s, e in full_exons if e == max(e2 for s2, e2 in full_exons)) + 1
                else:
                    junc = max(e for s, e in full_exons if s == min(s2 for s2, e2 in full_exons))
            row["exon_junc"] = str(junc)

            # Check if sc_alt or sc_ref equals exon_junc
            try:
                sc_alt_val = int(row["sc_alt"])
                if sc_alt_val == junc:
                    error_reason = f"sc_alt equals exon_junc ({junc})"
                    row.update({f: row.get(f, "") for f in main_fields if f not in row})
                    row["error_reason"] = error_reason
                    error_results.append(row)
                    continue
            except ValueError:
                pass

            try:
                sc_ref_val = int(row["sc_ref"])
                if sc_ref_val == junc:
                    error_reason = f"sc_ref equals exon_junc ({junc})"
                    row.update({f: row.get(f, "") for f in main_fields if f not in row})
                    row["error_reason"] = error_reason
                    error_results.append(row)
                    continue
            except ValueError:
                pass

           
            row["distance_sc_alt"] = str(calculate_transcript_distance(row["sc_alt"], alt_exons, strand, junc, full_exons, chrom, genome))

           
            row["distance_sc_ref"] = str(calculate_transcript_distance(row["sc_ref"], cds_exons, strand, junc, full_exons, chrom, genome))

           
            try:
                sc_alt_val = int(row["sc_alt"])
                sc_ref_val = int(row["sc_ref"])
                distance_sc_alt = int(row["distance_sc_alt"])
                distance_sc_ref = int(row["distance_sc_ref"])
                if distance_sc_alt == 0 and sc_alt_val != junc:
                    error_reason = f"distance_sc_alt is zero, but sc_alt ({sc_alt_val}) does not equal exon_junc ({junc})"
                    row.update({f: row.get(f, "") for f in main_fields if f not in row})
                    row["error_reason"] = error_reason
                    error_results.append(row)
                    continue
                if distance_sc_ref == 0 and sc_ref_val != junc:
                    error_reason = f"distance_sc_ref is zero, but sc_ref ({sc_ref_val}) does not equal exon_junc ({junc})"
                    row.update({f: row.get(f, "") for f in main_fields if f not in row})
                    row["error_reason"] = error_reason
                    error_results.append(row)
                    continue
            except ValueError:
                pass

           
            try:
                sc_alt_val = int(row["sc_alt"])
                sc_ref_val = int(row["sc_ref"])
                if row["Group"] == "Triplets" and row["In_frame_sc"] == "NO" and sc_alt_val != sc_ref_val:
                    error_reason = f"Triplets group has different sc_alt ({sc_alt_val}) and sc_ref ({sc_ref_val}) with In_frame_sc NO"
                    row.update({f: row.get(f, "") for f in main_fields if f not in row})
                    row["error_reason"] = error_reason
                    error_results.append(row)
                    continue
                if row["Group"] == "Non-Triplets" and sc_alt_val == sc_ref_val:
                    error_reason = f"Non-Triplets group has same sc_alt and sc_ref ({sc_alt_val})"
                    row.update({f: row.get(f, "") for f in main_fields if f not in row})
                    row["error_reason"] = error_reason
                    error_results.append(row)
                    continue
            except ValueError:
                pass

            
            end_UTR = 0
            try:
                if strand == "+":
                    end_UTR = max(e for s, e in full_exons)
                else:
                    end_UTR = min(s for s, e in full_exons)
                row["end_UTR"] = str(end_UTR)
            except ValueError:
                end_UTR = 0
                row["end_UTR"] = ""

            
            row["dist_ref_UTR"] = str(calculate_utr_distance(row["sc_ref"], cds_exons, strand, end_UTR, full_exons, chrom, genome))

            row["dist_alt_UTR"] = str(calculate_utr_distance(row["sc_alt"], alt_exons, strand, end_UTR, full_exons, chrom, genome))

            
            w_main.writerow(row)

    
    try:
        with open(false_exons_tsv, "w", newline="") as fex_out:
            w_false_exons = csv.DictWriter(fex_out, fieldnames=false_exons_fields, delimiter="\t")
            w_false_exons.writeheader()
            for excluded in false_exons_results:
                w_false_exons.writerow(excluded)
    except Exception as e:
        raise Exception(f"Error writing false_exons.tsv: {str(e)}")

   
    try:
        with open(error_tsv, "a", newline="") as error_out:
            w_error = csv.DictWriter(error_out, fieldnames=error_fields, delimiter="\t")
            for error_row in error_results:
                w_error.writerow(error_row)
    except Exception as e:
        raise Exception(f"Error writing error.tsv: {str(e)}")

except Exception as e:
    raise Exception(f"Error writing output files: {str(e)}")

print(f"Done. Outputs: {tsv_output}, {outcds_tsv}, {false_exons_tsv}, {error_tsv}, plus FASTAs in {sim_seq_dir}/")

