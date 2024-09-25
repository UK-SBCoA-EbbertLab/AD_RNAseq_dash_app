import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Import functions from our plotly_ggtranscript package
import sys
import os
sys.path.append("C:/Users/local_bag222/Desktop/dash_apps/AD_RNAseq_dash_app")
from plotly_ggtranscript import geom_range, geom_intron, to_intron, shorten_gaps, set_axis, calculate_cds_exon_difference, rescale_cds

# Create the sod1_annotation DataFrame
# Define the corrected data
data = {
    'seqnames': ['21'] * 35,
    'start': [
        31659666, 31659666, 31659770, 31663790, 31663790, 31666449, 31666449,
        31667258, 31667258, 31668471, 31668471, 31659693, 31659693, 31661549,
        31663790, 31666449, 31667258, 31659693, 31659693, 31659770, 31663790,
        31663790, 31666449, 31666449, 31667258, 31667258, 31668471, 31668471,
        31659709, 31659709, 31663790, 31666449, 31667258, 31668471
    ],
    'end': [
        31668931, 31659784, 31659784, 31663886, 31663886, 31666518, 31666518,
        31667375, 31667375, 31668931, 31668575, 31667341, 31659841, 31661734,
        31663886, 31666518, 31667341, 31668931, 31659841, 31659841, 31663886,
        31663886, 31666518, 31666518, 31667375, 31667375, 31668931, 31668575,
        31668931, 31660708, 31663886, 31666518, 31667375, 31668931
    ],
    'strand': ['+'] * 35,
    'type': [
        'transcript', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS', 'exon',
        'CDS', 'exon', 'CDS', 'transcript', 'exon', 'exon', 'exon', 'exon',
        'exon', 'transcript', 'exon', 'CDS', 'exon', 'CDS', 'exon', 'CDS',
        'exon', 'CDS', 'exon', 'CDS', 'transcript', 'exon', 'exon', 'exon',
        'exon', 'exon'
    ],
    'gene_name': ['SOD1'] * 35,
    'transcript_name': [
        'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202',
        'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202',
        'SOD1-202', 'SOD1-204', 'SOD1-204', 'SOD1-204', 'SOD1-204',
        'SOD1-204', 'SOD1-204', 'SOD1-201', 'SOD1-201', 'SOD1-201',
        'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-201',
        'SOD1-201', 'SOD1-201', 'SOD1-201', 'SOD1-203', 'SOD1-203',
        'SOD1-203', 'SOD1-203', 'SOD1-203', 'SOD1-203'
    ],
    'transcript_biotype': [
        'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding',
        'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding',
        'protein_coding', 'protein_coding', 'protein_coding', 'processed_transcript',
        'processed_transcript', 'processed_transcript', 'processed_transcript',
        'processed_transcript', 'processed_transcript', 'protein_coding',
        'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding',
        'protein_coding', 'protein_coding', 'protein_coding', 'protein_coding',
        'protein_coding', 'processed_transcript', 'processed_transcript',
        'processed_transcript', 'processed_transcript', 'processed_transcript',
        'processed_transcript'
    ],
    'exon_number': [
        np.nan, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0,
        np.nan, 1.0, 2.0, 3.0, 4.0, 5.0, np.nan, 1.0, 1.0, 2.0, 2.0,
        3.0, 3.0, 4.0, 4.0, 5.0, 5.0, np.nan, 1.0, 2.0, 3.0, 4.0,
        5.0
    ]
}

# Create the DataFrame
sod1_annotation = pd.DataFrame(data)

# Define a mapping from transcript_biotype to colors
biotype_colors = {
    'protein_coding': 'blue',
    'processed_transcript': 'green',
    'NA': 'gray'
}

## Add biotype colors
sod1_annotation['fillcolor'] = sod1_annotation['transcript_biotype'].map(biotype_colors)

# Extract exons
sod1_exons = sod1_annotation[sod1_annotation['type'] == 'exon']

## Extract CDS
sod1_cds = sod1_annotation[sod1_annotation['type'] == 'CDS']

## Calculate CDS and exon differences
sod1_cds_diff = calculate_cds_exon_difference(sod1_cds, sod1_exons)


# Use .loc to avoid SettingWithCopyWarning
sod1_exons = sod1_exons.copy()

## Rescale SOD exons
sod1_rescaled = shorten_gaps(exons=sod1_exons, introns=to_intron(sod1_exons, "transcript_name"), group_var = "transcript_name")

## Define rescaled exons and introns
sod1_rescaled_exons = sod1_rescaled.loc[sod1_rescaled["type"] == "exon"].copy()
sod1_rescaled_introns = sod1_rescaled.loc[sod1_rescaled["type"] == "intron"].copy()

## Correct CDS coordinates for new exon coordinates
sod1_rescaled_cds = rescale_cds(sod1_cds_diff, sod1_rescaled_exons)

# Create the plot
# fig = go.Figure()

## Debug
print(sod1_cds_diff.head())
print(sod1_rescaled_cds.loc[sod1_rescaled_cds["transcript_name"] == "SOD1-202"].head(10))
print(sod1_rescaled_exons.loc[sod1_rescaled_exons["transcript_name"] == "SOD1-202"].head(10))

# # Add exons using geom_range, passing the fillcolor directly
# exon_traces = geom_range(
#     data=sod1_rescaled_exons,
#     x_start='start',
#     x_end='end',
#     y='transcript_name',
#     fill=sod1_rescaled_exons['fillcolor'], 
#     height=0.5
# )

# ## Add CDS traces
# cds_traces = geom_range(
#     data=sod1_rescaled_cds,
#     x_start='start',
#     x_end='end',
#     y='transcript_name',
#     fill=sod1_rescaled_cds['fillcolor'],
#     height= 0.25
# )


# # Create introns and add them using geom_intron
# #sod1_introns = to_intron(sod1_exons, group_var="transcript_name")
# intron_traces = geom_intron(
#     data=sod1_rescaled_introns,
#     x_start='start',
#     x_end='end',
#     y='transcript_name',
#     strand='strand'
# )

# # Add exons, CDS, and introns as before
# for trace in exon_traces:
#     fig.add_shape(trace)

# for trace in cds_traces:
#     fig.add_shape(trace)

# for trace in intron_traces:
#     if isinstance(trace, dict):
#         fig.add_shape(trace)
#     else:
#         fig.add_trace(trace)

# # Call the new function to set the genomic axis range
# fig = set_axis(fig, sod1_rescaled_exons, sod1_rescaled_introns)

# # Update layout and show the plot
# fig.update_layout(
#     title="SOD1 Transcript Structure",
#     xaxis_title="Genomic Position",
#     yaxis_title="Transcript",
#     height=400,
#     width=800,
#     showlegend=False
# )

# # Show or save the plot
# fig.show()
# fig.write_html("sod1_transcript_structure.html")
