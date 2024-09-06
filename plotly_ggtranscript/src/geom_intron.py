import plotly.graph_objects as go
import pandas as pd
import numpy as np

def geom_intron(data, x_start='start', x_end='end', y='transcript_name', 
                strand='strand', color='black', line_width=0.5, opacity=1,
                arrow_color='black', arrow_size=0.3, arrow_min_intron_length=0):
    """
    Create a Plotly figure with introns represented as lines with strand arrows.

    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame containing the intron data.
    x_start : str, optional
        Column name for the start position of introns. Default is 'start'.
    x_end : str, optional
        Column name for the end position of introns. Default is 'end'.
    y : str, optional
        Column name for the y-axis values (e.g., transcript names). Default is 'transcript_name'.
    strand : str, optional
        Column name for the strand information. Default is 'strand'.
    color : str, optional
        Color of the intron lines. Default is 'black'.
    line_width : float, optional
        Width of the intron lines. Default is 0.5.
    opacity : float, optional
        Opacity of the intron lines. Default is 1.
    arrow_color : str, optional
        Color of the strand arrows. Default is 'black'.
    arrow_size : float, optional
        Size of the strand arrows. Default is 0.3.
    arrow_min_intron_length : int, optional
        Minimum intron length to draw an arrow. Default is 0.

    Returns:
    --------
    plotly.graph_objects.Figure
        A Plotly figure object containing the intron plot with strand arrows.
    """
    fig = go.Figure()

    # Calculate y positions
    unique_y = data[y].unique()
    y_dict = {val: i for i, val in enumerate(unique_y)}

    # Add intron lines
    for idx, row in data.iterrows():
        y_pos = y_dict[row[y]]
        fig.add_shape(
            type="line",
            x0=row[x_start],
            x1=row[x_end],
            y0=y_pos,
            y1=y_pos,
            line=dict(color=color, width=line_width),
            opacity=opacity
        )

        # Add strand arrows
        intron_length = abs(row[x_end] - row[x_start])
        if intron_length > arrow_min_intron_length:
            mid_point = (row[x_start] + row[x_end]) / 2
            arrow_direction = 1 if row[strand] == '+' else -1
            fig.add_annotation(
                x=mid_point,
                y=y_pos,
                ax=mid_point + arrow_direction * intron_length * 0.05,
                ay=y_pos,
                xref="x",
                yref="y",
                axref="x",
                ayref="y",
                showarrow=True,
                arrowhead=2,
                arrowsize=arrow_size,
                arrowwidth=line_width,
                arrowcolor=arrow_color
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
    #    'strand': ['+', '+', '+', '-', '-']
    #})

    # Create the plot
    #fig = geom_intron(data, arrow_min_intron_length=30)
    #fig.show()