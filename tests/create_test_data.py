import pandas as pd

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

df = pd.DataFrame(data)
df = df.loc[df["type"].isin(["gene", "transcript", "exon", "CDS"])].copy()
df['exon_number'] = pd.NA

transcripts = df['transcript_name'].unique()
transcripts = [t for t in transcripts if t != 'NA']

for transcript in transcripts:
    # Get indices of this transcript
    indices_transcript = df[df['transcript_name'] == transcript].index
    
    # Get exons
    df_exons = df.loc[indices_transcript]
    df_exons = df_exons[df_exons['type'] == 'exon']
    # Sort exons by start position
    df_exons = df_exons.sort_values(by='start')
    # Assign exon numbers starting from 1
    exon_numbers = list(range(1, len(df_exons) + 1))
    df_exons = df_exons.assign(exon_number=exon_numbers)
    # Update exon_number in the main dataframe for exons
    df.loc[df_exons.index, 'exon_number'] = exon_numbers
    
    # Map from exon start/end to exon number
    exon_coords = list(zip(df_exons['start'], df_exons['end'], exon_numbers))
    
    # Now, for each row in df_transcript, assign exon_number if applicable
    for index in indices_transcript:
        row = df.loc[index]
        if row['type'] in ['CDS', 'start_codon', 'stop_codon', 'five_prime_utr', 'three_prime_utr']:
            # Find the exon where the row's start position is within exon start and end
            assigned = False
            for exon_start, exon_end, exon_num in exon_coords:
                if exon_start <= row['start'] <= exon_end:
                    df.loc[index, 'exon_number'] = exon_num
                    assigned = True
                    break
            if not assigned:
                df.loc[index, 'exon_number'] = pd.NA
        elif row['type'] == 'transcript':
            # Set exon_number to NA
            df.loc[index, 'exon_number'] = pd.NA

# Display the final DataFrame
print(df)

## Output test data
df.to_csv("../raw_data/test_data/test_annotations.tsv", index=False, sep="\t")