import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Import functions from our plotly_ggtranscript package
import sys
import os
sys.path.append("C:/Users/local_bag222/Desktop/dash_apps/AD_RNAseq_dash_app")
from plotly_ggtranscript import geom_range, geom_intron, to_intron

# Create the sod1_annotation DataFrame
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

sod1_annotation = pd.DataFrame(data)
sod1_annotation.dropna(inplace=True)

# Define a mapping from transcript_biotype to colors
biotype_colors = {
    'protein_coding': 'blue',
    'processed_transcript': 'green',
    'NA': 'gray'
}

# Extract exons
sod1_exons = sod1_annotation[sod1_annotation['type'] == 'exon']

# Use .loc to avoid SettingWithCopyWarning
sod1_exons.loc[:, 'fillcolor'] = sod1_exons['transcript_biotype'].map(biotype_colors)

# Create the plot
fig = go.Figure()

# Add exons using geom_range, passing the fillcolor directly
exon_traces = geom_range(
    data=sod1_exons,
    x_start='start',
    x_end='end',
    y='transcript_name',
    fill=sod1_exons['fillcolor']  # Pass the color-mapped column
)
for trace in exon_traces:
    fig.add_shape(trace)  # Add shapes to the figure using add_shape

# Create introns and add them using geom_intron
sod1_introns = to_intron(sod1_exons, group_var="transcript_name")
intron_traces = geom_intron(
    data=sod1_introns,
    x_start='start',
    x_end='end',
    y='transcript_name',
    strand='strand'
)

# Add intron shapes and arrows
for trace in intron_traces:
    if isinstance(trace, dict):  # Add shapes
        fig.add_shape(trace)
    else:  # Add scatter traces
        fig.add_trace(trace)

# Update layout
fig.update_layout(
    title="SOD1 Transcript Structure",
    xaxis_title="Genomic Position",
    yaxis_title="Transcript",
    height=400,
    width=800,
    showlegend=False
)

# Show the plot
fig.show()

# If you want to save the plot as an HTML file
fig.write_html("sod1_transcript_structure.html")

print("Plot has been generated and saved as 'sod1_transcript_structure.html'")