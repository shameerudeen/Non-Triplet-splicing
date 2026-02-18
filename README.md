# Non-Triplet-splicing
# CDS-Level Alternative Splicing Simulation & NMD Annotation Pipeline (see psuedocode in the methods section for additional explanation)

This pipeline consists of three sequential Python scripts that:

1. Assign a reference transcript per gene (rt.py)
2. Determine whether a splicing event is included or skipped within the CDS (status.py)
3. Simulate alternative splicing at the CDS level and compute NMD-related features (simulation.py)

The workflow is designed for transcriptome-wide alternative splicing analysis downstream of event detection tools (e.g., JUM, MAJIQ and rMATS).


SYSTEM REQUIREMENTS

Operating Systems Tested:
- CentOS 7 (HPC environment) (any system with below mentioned tools can be used)

Tested Software Environment:
- Python 3.9.18 (conda-forge)
- GCC 12.3.0
- Biopython 1.85
- pandas ≥ 1.5

Genome Annotation Used:
- WormBase WS290 (#!genebuild-version WS290)

Required Python Packages:
- biopython==1.85
- pandas

Install with:

conda create -n splicing_env python=3.9.18 biopython=1.85 pandas
conda activate splicing_env

Typical installation time:
< 2 minutes


PIPELINE OVERVIEW

Input:
- JUM splicing event TSV
- WormBase WS290 annotation.gtf
- Matching genome.fa

Step 1 → rt.py
Step 2 → status.py
Step 3 → simulation.py


STEP 1: ASSIGN REFERENCE TRANSCRIPT (rt.py)

Purpose:
For each gene, identify the transcript with the longest CDS span.
This transcript is assigned as the reference_transcript.

Input:
- TSV containing:
  - Gene
  - start
  - end
- annotation.gtf (WS290)

Output:
- TSV identical to input, plus:
  - reference_transcript column

Run:

python rt.py


STEP 2: DETERMINE INCLUSION STATUS (status.py)

Purpose:
Determine whether the splicing event coordinates fall entirely within
a CDS exon of the reference transcript.

Input:
- Output from rt.py
- annotation.gtf

Output:
- TSV with added:
  - status column ("included" or "skipped")

Logic:
- Extract CDS exons for reference transcript.
- If event start–end is fully within a CDS exon → "included"
- Otherwise → "skipped"

Run:

python status.py


STEP 3: CDS SIMULATION & NMD ANNOTATION (simulation.py)

Purpose:
Simulate alternative exon inclusion/exclusion and compute:

- sc_ref (reference stop codon)
- sc_alt (alternate stop codon)
- In_frame_sc
- recoded_exons
- recoded nucleotide length
- dual_coding length
- distance to exon junction
- distance to 3' UTR
- downstream exon counts

Also reconstruct:
- Full reference CDS sequence
- Full alternate CDS sequence
- Simulated region sequence

Input:
- Output from status.py
- annotation.gtf (WS290)
- genome.fa (matching WS290)

Outputs:

1. Main annotated TSV
2. Events outside CDS boundaries
3. Events overlapping multiple CDS exons
4. Logical inconsistencies (error file)
5. FASTA directory containing:
   - simulated region sequence
   - region length
   - alternate CDS
   - reference CDS

Run:

python simulation.py

Expected runtime:
~100 events → < 10 seconds
~5,000 events → < 5 minutes (C. elegans)


INPUT REQUIREMENTS

Required TSV columns (initial JUM output):
- Gene
- start
- end
- chr
- strand
  
Additional columns retained if present:

- avgwtPSI
- avgsmgPSI
- deltaPSI
- Group
- ID
- category
- Phase
- AS_event_ID
- NMD
- Significance

COLUMNS ADDED DURING PIPELINE ( see psuedocode in the methods section for additional explanation regarding the meaning of each col)

Step 1 (rt.py):
- reference_transcript

Step 2 (status.py):
- status

Step 3 (simulation.py):

Stop codon and frame:
- sc_ref
- sc_alt
- In_frame_sc

Exon structure:
- exon_after_end
- exon_after_sc_alt
- exon_junc

Recoding metrics:
- recoded_exons
- ntd_recoded_ref
- ntds_altv_rec
- dual_coding

Distance metrics:
- distance_sc_ref
- distance_sc_alt
- end_UTR
- dist_ref_UTR
- dist_alt_UTR

These columns are appended to the input table and written to the main annotated TSV.


REPRODUCIBILITY

To reproduce results exactly:

- Use WormBase WS290 annotation
- Use matching WS290 genome FASTA
- Use Python 3.9.18
- Use Biopython 1.85
- Use pandas ≥ 1.5


ERROR HANDLING

Events are excluded if:

- Outside CDS boundaries
- No CDS annotation
- Overlapping multiple CDS exons
- Stop codon equals exon junction
- Logical inconsistency

Excluded events are written to separate TSV files.

