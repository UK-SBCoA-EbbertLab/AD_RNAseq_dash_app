import plotly.graph_objects as go
import pandas as pd
import numpy as np

def geom_range(data, x_start='start', x_end='end', y='transcript_name',
               fill='grey', color='black', opacity=1, line_width=0.25, height=0.5):
    """
    Create a list of Plotly range (exon) trace objects.
    
    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame containing the range (exon) data.
    x_start : str, optional
        Column name for the start position of exons. Default is 'start'.
    x_end : str, optional
        Column name for the end position of exons. Default is 'end'.
    y : str, optional
        Column name for the y-axis values (e.g., transcript names). Default is 'transcript_name'.
    fill : str or list, optional
        Fill color or a list of colors. Default is 'grey'.
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
    list
        A list of Plotly trace objects representing the exons.
    """
    traces = []
    unique_y = data[y].unique()
    y_dict = {val: i for i, val in enumerate(unique_y)}
    
    if isinstance(fill, str):
        fill = [fill] * len(data)
    
    for idx, row in data.iterrows():
        y_pos = y_dict[row[y]]
        trace = dict(
            type="rect",
            x0=row[x_start],
            x1=row[x_end],
            y0=y_pos - height/2,
            y1=y_pos + height/2,
            fillcolor=fill[idx],
            line=dict(color=color, width=line_width),
            opacity=opacity
        )
        traces.append(trace)
    
    return traces
