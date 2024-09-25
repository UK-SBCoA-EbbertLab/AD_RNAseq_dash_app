## Import libraries and set path
import pandas as pd
import os
import sys
sys.path.append("C:/Users/local_bag222/Desktop/dash_apps/AD_RNAseq_dash_app")

## Define data
data = {
    'seqnames': ['21'] * 43,
    'start': [31659666, 31659666, 31659666, 31659770, 31659770, 31663790, 31663790, 31666449, 31666449, 31667258, 31667258, 31668471, 31668471, 31668576, 31659666, 31668579, 31659693, 31659693, 31661549, 31663790, 31666449, 31667258, 31659693, 31659693, 31659770, 31659770, 31663790, 31663790, 31666449, 31666449, 31667258, 31667258, 31668471, 31668471, 31668576, 31659693, 31668579, 31659709, 31659709, 31663790, 31666449, 31667258, 31668471],
    'end': [31668931, 31668931, 31659784, 31659784, 31659772, 31663886, 31663886, 31666518, 31666518, 31667375, 31667375, 31668931, 31668575, 31668578, 31659769, 31668931, 31667341, 31659841, 31661734, 31663886, 31666518, 31667341, 31668931, 31659841, 31659841, 31659772, 31663886, 31663886, 31666518, 31666518, 31667375, 31667375, 31668931, 31668575, 31668578, 31659769, 31668931, 31668931, 31660708, 31663886, 31666518, 31667375, 31668931],
    'strand': ['+'] * 43,
    'type': ['gene', 'transcript', 'exon', 'CDS', 'start_codon', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS', 'stop_codon', 'five_prime_utr', 'three_prime_utr', 'transcript', 'exon', 'exon', 'exon', 'exon', 'exon', 'transcript', 'exon', 'CDS', 'start_codon', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS', 'stop_codon', 'five_prime_utr', 'three_prime_utr', 'transcript', 'exon', 'exon', 'exon', 'exon', 'exon'],
    'gene_name': ['SOD1'] * 43,
    'transcript_name': ['NA', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-204', 'SOD1-204', 'SOD1-204', 'SOD1-204', 'SOD1-204', 'SOD1-204', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-203', 'SOD1-203', 'SOD1-203', 'SOD1-203', 'SOD1-203', 'SOD1-203'],
    'transcript_biotype': ['NA', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript', 'processed_transcript']
}

# Create a DataFrame from the input data
df = pd.DataFrame(data)

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
    
    # Extract exon rows for the current transcript
    df_exons = df.loc[indices_transcript]
    df_exons = df_exons[df_exons['type'] == 'exon']
    
    # Sort exons by their start position (assuming exon numbers are ordered by genomic location)
    df_exons = df_exons.sort_values(by='start')
    
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
        
        # For specific row types (CDS, start/stop codon, UTRs), try to assign an exon number
        if row['type'] in ['CDS', 'start_codon', 'stop_codon', 'five_prime_utr', 'three_prime_utr']:
            assigned = False
            # Assign the exon number if the row's 'start' position falls within an exon range
            for exon_start, exon_end, exon_num in exon_coords:
                if exon_start <= row['start'] <= exon_end:
                    df.loc[index, 'exon_number'] = exon_num
                    assigned = True
                    break
            # If no exon matches the row's position, set 'exon_number' to NA
            if not assigned:
                df.loc[index, 'exon_number'] = pd.NA
        
        # For transcript rows, we explicitly set 'exon_number' to NA (not applicable)
        elif row['type'] == 'transcript':
            df.loc[index, 'exon_number'] = pd.NA

# Output the updated DataFrame to a test file in TSV format
df.to_csv("./raw_data/test_data/test_annotations.tsv", index=False, sep="\t")