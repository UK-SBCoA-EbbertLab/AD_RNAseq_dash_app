import plotly.graph_objects as go
import pandas as pd
import numpy as np

def geom_range(data, x_start='start', x_end='end', y='transcript_name', 
               fill='grey', color='black', opacity=1, line_width=0.25, 
               height=0.5):
    """
    Create a Plotly figure with genomic ranges represented as rectangles.

    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame containing the genomic range data.
    x_start : str, optional
        Column name for the start position of ranges. Default is 'start'.
    x_end : str, optional
        Column name for the end position of ranges. Default is 'end'.
    y : str, optional
        Column name for the y-axis values (e.g., transcript names). Default is 'transcript_name'.
    fill : str or list, optional
        Fill color for the rectangles. Can be a single color or a list of colors. Default is 'grey'.
    color : str, optional
        Color of the rectangle borders. Default is 'black'.
    opacity : float, optional
        Opacity of the rectangles. Default is 1.
    line_width : float, optional
        Width of the rectangle borders. Default is 0.25.
    height : float, optional
        Height of the rectangles. Default is 0.5.

    Returns:
    --------
    plotly.graph_objects.Figure
        A Plotly figure object containing the genomic range plot.
    """
    fig = go.Figure()

    # Calculate y positions
    unique_y = data[y].unique()
    y_dict = {val: i for i, val in enumerate(unique_y)}

    # If fill is a single color, repeat it for all ranges
    if isinstance(fill, str):
        fill = [fill] * len(data)

    for idx, row in data.iterrows():
        y_pos = y_dict[row[y]]
        fig.add_shape(
            type="rect",
            x0=row[x_start],
            x1=row[x_end],
            y0=y_pos - height/2,
            y1=y_pos + height/2,
            fillcolor=fill[idx],
            line=dict(color=color, width=line_width),
            opacity=opacity
        )

    # Set layout
    fig.update_layout(
        yaxis=dict(
            tickmode='array',
            tickvals=list(y_dict.values()),
            ticktext=list(y_dict.keys()),
            zeroline=False,
        ),
        xaxis=dict(title="Genomic Position"),
        showlegend=False
    )

    return fig

# Example usage:
#if __name__ == "__main__":
    # Create sample data
    #data = pd.DataFrame({
    #    'start': [100, 200, 300, 150, 250],
    #    'end': [150, 250, 350, 200, 300],
    #    'transcript_name': ['Gene1', 'Gene1', 'Gene1', 'Gene2', 'Gene2'],
    #    'feature': ['exon1', 'exon2', 'exon3', 'exon1', 'exon2']
    #})

    # Create the plot
    #fig = geom_range(data, fill=data['feature'].map({'exon1': 'red', 'exon2': 'green', 'exon3': 'blue'}))
    #fig.show()