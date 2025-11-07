import pandas as pd
from collections import defaultdict

def parse_gtf(gtf_file):
    """Parse GTF file to extract first exon information for each transcript."""
    gene_transcripts = defaultdict(list)
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] == 'exon':
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]
                # Extract gene_id and transcript_id
                attr_dict = {}
                for attr in attributes.split(';'):
                    if attr.strip():
                        key, value = attr.strip().split(' ', 1)
                        attr_dict[key] = value.strip('"')
                gene_id = attr_dict.get('gene_id')
                transcript_id = attr_dict.get('transcript_id')
                exon_number = attr_dict.get('exon_number')
                if gene_id and transcript_id and exon_number == '1':
                    gene_transcripts[gene_id].append({
                        'transcript_id': transcript_id,
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand
                    })
    return gene_transcripts

def get_transcripts_in_first_exon(row, gene_transcripts):
    """Get transcripts where region overlaps or spans the first exon."""
    gene_id = row['ID']
    region_start = int(row['start'])
    region_end = int(row['end'])
    if gene_id not in gene_transcripts:
        return None
    transcripts = []
    for transcript in gene_transcripts[gene_id]:
        exon_start = transcript['start']
        exon_end = transcript['end']
        if (region_start >= exon_start and region_end <= exon_end) or \
           (region_start <= exon_start and region_end >= exon_end):
            transcripts.append(transcript['transcript_id'])
    return ','.join(transcripts) if transcripts else None

def main(casst_file, gtf_file, output_file):
    casst_df = pd.read_csv(casst_file, sep='\t')
    gene_transcripts = parse_gtf(gtf_file)
    casst_df['Transcripts'] = casst_df.apply(lambda row: get_transcripts_in_first_exon(row, gene_transcripts), axis=1)
    output_df = casst_df[casst_df['Transcripts'].notnull()]
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    casst_file = 'tsv containing JUM output coordinates'
    gtf_file = 'annotation.gtf'
    output_file = 'non_triplets_first_exon.txt'
    main(casst_file, gtf_file, output_file)
