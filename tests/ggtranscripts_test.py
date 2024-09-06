import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Import functions from our plotly_ggtranscript package
from plotly_ggtranscript.src import geom_range, geom_intron, to_intron

# Load example data (you'll need to replace this with actual data loading)
# For this example, we'll create a sample DataFrame similar to sod1_annotation
sod1_annotation = pd.DataFrame({
    'seqnames': ['21'] * 6,
    'start': [31659666, 31659666, 31659666, 31659770, 31659770, 31663790],
    'end': [31668931, 31668931, 31659784, 31659784, 31659772, 31663886],
    'strand': ['+'] * 6,
    'type': ['gene', 'transcript', 'exon', 'CDS', 'start_codon', 'exon'],
    'gene_name': ['SOD1'] * 6,
    'transcript_name': [None, 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202', 'SOD1-202'],
    'transcript_biotype': ['protein_coding'] * 6
})

print("First few rows of sod1_annotation:")
print(sod1_annotation.head())

# Extract exons
sod1_exons = sod1_annotation[sod1_annotation['type'] == 'exon']

# Create introns
sod1_introns = to_intron(sod1_exons, group_var="transcript_name")

# Create the plot
fig = make_subplots(rows=1, cols=1, subplot_titles=['SOD1 Transcript Structure'])

# Add exons using geom_range
exon_traces = geom_range(
    data=sod1_exons,
    x_start='start',
    x_end='end',
    y='transcript_name',
    fill=sod1_exons['transcript_biotype']
)
for trace in exon_traces:
    fig.add_trace(trace)

# Add introns using geom_intron
intron_traces = geom_intron(
    data=sod1_introns,
    x_start='start',
    x_end='end',
    y='transcript_name',
    strand='strand'
)
for trace in intron_traces:
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