## Import libraries and set path
import pandas as pd
import os
import sys
sys.path.append("C:/Users/local_bag222/Desktop/dash_apps/AD_RNAseq_dash_app")

## Define data
sod1_data = {
    'seqnames': ['21'] * 43,
    'start': [31659666, 31659666, 31659666, 31659770, 31659770, 31663790, 31663790, 31666449, 31666449, 31667258, 31667258, 31668471, 31668471, 31668576, 31659666, 31668579, 31659693, 31659693, 31661549, 31663790, 31666449, 31667258, 31659693, 31659693, 31659770, 31659770, 31663790, 31663790, 31666449, 31666449, 31667258, 31667258, 31668471, 31668471, 31668576, 31659693, 31668579, 31659709, 31659709, 31663790, 31666449, 31667258, 31668471],
    'end': [31668931, 31668931, 31659784, 31659784, 31659772, 31663886, 31663886, 31666518, 31666518, 31667375, 31667375, 31668931, 31668575, 31668578, 31659769, 31668931, 31667341, 31659841, 31661734, 31663886, 31666518, 31667341, 31668931, 31659841, 31659841, 31659772, 31663886, 31663886, 31666518, 31666518, 31667375, 31667375, 31668931, 31668575, 31668578, 31659769, 31668931, 31668931, 31660708, 31663886, 31666518, 31667375, 31668931],
    'strand': ['+'] * 43,
    'type': ['gene', 'transcript', 'exon', 'CDS', 'start_codon', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS', 'stop_codon', 'five_prime_utr', 'three_prime_utr', 'transcript', 'exon', 'exon', 'exon', 'exon', 'exon', 'transcript', 'exon', 'CDS', 'start_codon', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS', 'stop_codon', 'five_prime_utr', 'three_prime_utr', 'transcript', 'exon', 'exon', 'exon', 'exon', 'exon'],
    'gene_name': ['SOD1'] * 43,
    'transcript_name': ['NA', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-204', 'SOD1-204', 'SOD1-204', 'SOD1-204', 'SOD1-204', 'SOD1-204', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-203', 'SOD1-203', 'SOD1-203', 'SOD1-203', 'SOD1-203', 'SOD1-203'],
    'transcript_biotype': ['NA', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript']
}

# New data for RUNX1 gene (negative strand)
runx1_data = {
    'seqnames': ['21'] * 16,
    'start': [
        34620210,  # gene
        34620210,  # transcript RUNX1-201
        34845000,  # exon 1
        34750000,  # exon 2
        34660000,  # exon 3
        34620210,  # exon 4
        34845000,  # CDS 1
        34750000,  # CDS 2
        34660000,  # CDS 3
        34660000,  # start_codon
        34845523,  # stop_codon
        34845526,  # five_prime_utr
        34659500,  # three_prime_utr
        34620210,  # transcript RUNX1-202
        34800000,  # exon 1
        34620210   # exon 2
    ],
    'end': [
        34845525,  # gene
        34845525,  # transcript RUNX1-201
        34845525,  # exon 1
        34750500,  # exon 2
        34660500,  # exon 3
        34620500,  # exon 4
        34845525,  # CDS 1
        34750500,  # CDS 2
        34660500,  # CDS 3
        34660002,  # start_codon
        34845525,  # stop_codon
        34846000,  # five_prime_utr
        34660000,  # three_prime_utr
        34800000,  # transcript RUNX1-202
        34800500,  # exon 1
        34620500   # exon 2
    ],
    'strand': ['-'] * 16,
    'type': [
        'gene',            # gene
        'transcript',      # transcript RUNX1-201
        'exon',            # exon 1
        'exon',            # exon 2
        'exon',            # exon 3
        'exon',            # exon 4
        'CDS',             # CDS 1
        'CDS',             # CDS 2
        'CDS',             # CDS 3
        'start_codon',     # start_codon
        'stop_codon',      # stop_codon
        'five_prime_utr',  # five_prime_utr
        'three_prime_utr', # three_prime_utr
        'transcript',      # transcript RUNX1-202
        'exon',            # exon 1
        'exon'             # exon 2
    ],
    'gene_name': ['RUNX1'] * 16,
    'transcript_name': [
        'NA',             # gene
        'RUNX1-201',      # transcript RUNX1-201
        'RUNX1-201',      # exon 1
        'RUNX1-201',      # exon 2
        'RUNX1-201',      # exon 3
        'RUNX1-201',      # exon 4
        'RUNX1-201',      # CDS 1
        'RUNX1-201',      # CDS 2
        'RUNX1-201',      # CDS 3
        'RUNX1-201',      # start_codon
        'RUNX1-201',      # stop_codon
        'RUNX1-201',      # five_prime_utr
        'RUNX1-201',      # three_prime_utr
        'RUNX1-202',      # transcript RUNX1-202
        'RUNX1-202',      # exon 1
        'RUNX1-202'       # exon 2
    ],
    'transcript_biotype': [
        'NA',                   # gene
        'protein_coding',       # transcript RUNX1-201
        'protein_coding',       # exon 1
        'protein_coding',       # exon 2
        'protein_coding',       # exon 3
        'protein_coding',       # exon 4
        'protein_coding',       # CDS 1
        'protein_coding',       # CDS 2
        'protein_coding',       # CDS 3
        'protein_coding',       # start_codon
        'protein_coding',       # stop_codon
        'protein_coding',       # five_prime_utr
        'protein_coding',       # three_prime_utr
        'processed_transcript', # transcript RUNX1-202
        'processed_transcript', # exon 1
        'processed_transcript'  # exon 2
    ]
}

# Create DataFrames from the input data
df_sod1 = pd.DataFrame(sod1_data)
df_runx1 = pd.DataFrame(runx1_data)

# Combine the DataFrames
df = pd.concat([df_sod1, df_runx1], ignore_index=True)

# Filter the DataFrame to keep only relevant rows (i.e., gene, transcript, exon, CDS types)
df = df.loc[df["type"].isin(["gene", "transcript", "exon", "CDS"])].copy()

# Initialize the 'exon_number' column as null (for future assignment)
df['exon_number'] = pd.NA

# Get the list of unique transcript names, excluding any entries labeled as 'NA'
transcripts = df['transcript_name'].unique()
transcripts = [t for t in transcripts if t != 'NA']

# Loop over each transcript to assign exon numbers
for transcript in transcripts:
    # Get the indices corresponding to the current transcript in the DataFrame
    indices_transcript = df[df['transcript_name'] == transcript].index

    # Get the strand for this transcript
    strand = df.loc[indices_transcript, 'strand'].unique()
    if len(strand) == 1:
        strand = strand[0]
    else:
        print(f"Warning: multiple strands found for transcript {transcript}")
        continue

    # Extract exon rows for the current transcript
    df_exons = df.loc[indices_transcript]
    df_exons = df_exons[df_exons['type'] == 'exon']

    # Sort exons by their start position, order depends on strand
    if strand == '+':
        df_exons = df_exons.sort_values(by='start', ascending=True)
    elif strand == '-':
        df_exons = df_exons.sort_values(by='start', ascending=False)
    else:
        print(f"Warning: unknown strand {strand} for transcript {transcript}")
        continue

    # Create a list of exon numbers starting from 1, matching the number of exons
    exon_numbers = list(range(1, len(df_exons) + 1))

    # Assign the calculated exon numbers to the 'exon_number' column of exon rows
    df_exons = df_exons.assign(exon_number=exon_numbers)

    # Update the main DataFrame with the assigned exon numbers for the exons
    df.loc[df_exons.index, 'exon_number'] = exon_numbers

    # Create a list of tuples containing the start, end, and exon number for the current transcript
    exon_coords = list(zip(df_exons['start'], df_exons['end'], exon_numbers))

    # Iterate over each row of the current transcript to assign exon numbers where applicable
    for index in indices_transcript:
        row = df.loc[index]

        # For specific row types, try to assign an exon number
        if row['type'] in ['CDS', 'start_codon', 'stop_codon', 'five_prime_utr', 'three_prime_utr']:
            assigned = False
            # Assign the exon number if the row's position overlaps with an exon range
            for exon_start, exon_end, exon_num in exon_coords:
                if row['start'] >= exon_start and row['end'] <= exon_end:
                    df.loc[index, 'exon_number'] = exon_num
                    assigned = True
                    break
            # If no exon matches the row's position, set 'exon_number' to NA
            if not assigned:
                df.loc[index, 'exon_number'] = pd.NA

        # For transcript and gene rows, set 'exon_number' to NA
        elif row['type'] in ['transcript', 'gene']:
            df.loc[index, 'exon_number'] = pd.NA

# Output the updated DataFrame to a test file in TSV format
df.to_csv("./raw_data/test_data/test_annotations.tsv", index=False, sep="\t")

# Optional: Print the DataFrame to verify
print(df)