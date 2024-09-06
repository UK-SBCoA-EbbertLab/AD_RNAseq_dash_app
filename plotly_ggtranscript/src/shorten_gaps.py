import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import List, Union

def shorten_gaps(exons: pd.DataFrame, introns: pd.DataFrame, 
                 group_var: Union[str, List[str]] = None, 
                 target_gap_width: int = 100) -> pd.DataFrame:
    """
    Improve transcript structure visualization by shortening gaps between exons.

    Parameters:
    -----------
    exons : pd.DataFrame
        DataFrame containing exon coordinates.
    introns : pd.DataFrame
        DataFrame containing intron coordinates.
    group_var : str or List[str], optional
        Column name(s) for grouping transcripts. Default is None.
    target_gap_width : int, optional
        The width in base pairs to shorten the gaps to. Default is 100.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the re-scaled coordinates of introns and exons
        with shortened gaps.
    """
    # Input checks
    for df in [exons, introns]:
        assert all(col in df.columns for col in ['start', 'end', 'strand', 'seqnames']), \
            "Both exons and introns DataFrames must have 'start', 'end', 'strand', and 'seqnames' columns"

    # Ensure 'type' column exists
    exons = _get_type(exons, "exons")
    introns = _get_type(introns, "introns")

    # Adjust intron boundaries
    introns = introns.copy()
    introns['start'] += 1
    introns['end'] -= 1

    # Get gaps (regions that don't overlap with any exons)
    gaps = _get_gaps(exons)

    # Map gaps to introns
    gap_map = _get_gap_map(introns, gaps)

    # Shorten gaps
    introns_shortened = _get_shortened_gaps(introns, gaps, gap_map, group_var, target_gap_width)

    # Handle transcript start gaps if multiple transcripts
    if group_var:
        tx_start_gaps = _get_tx_start_gaps(exons, group_var)
        gap_map_tx_start = _get_gap_map(tx_start_gaps, gaps)
        tx_start_gaps_shortened = _get_shortened_gaps(tx_start_gaps, gaps, gap_map_tx_start, group_var, target_gap_width)
        tx_start_gaps_shortened = tx_start_gaps_shortened.drop(columns=['start', 'end', 'strand', 'seqnames'])
    else:
        tx_start_gaps_shortened = None

    # Get rescaled transcripts
    rescaled_tx = _get_rescaled_txs(exons, introns_shortened, tx_start_gaps_shortened, group_var)

    return rescaled_tx

def _get_type(df: pd.DataFrame, df_type: str) -> pd.DataFrame:
    if 'type' not in df.columns:
        df = df.copy()
        df['type'] = 'exon' if df_type == 'exons' else 'intron'
    elif df_type == 'introns':
        assert all(df['type'] == 'intron'), "All 'type' values in introns DataFrame must be 'intron'"
    return df

def _get_gaps(exons: pd.DataFrame) -> pd.DataFrame:
    # Ensure single chromosome and strand
    assert exons['seqnames'].nunique() == 1 and exons['strand'].nunique() == 1, \
        "Exons must be from a single chromosome and strand"

    # Sort and merge overlapping exons
    exons_sorted = exons.sort_values('start')
    merged_exons = []
    for _, exon in exons_sorted.iterrows():
        if not merged_exons or exon['start'] > merged_exons[-1]['end']:
            merged_exons.append({'start': exon['start'], 'end': exon['end']})
        else:
            merged_exons[-1]['end'] = max(merged_exons[-1]['end'], exon['end'])

    # Find gaps between merged exons
    gaps = []
    for i in range(1, len(merged_exons)):
        gap_start = merged_exons[i-1]['end'] + 1
        gap_end = merged_exons[i]['start'] - 1
        if gap_start <= gap_end:
            gaps.append({'start': gap_start, 'end': gap_end})

    return pd.DataFrame(gaps)

def _get_tx_start_gaps(exons: pd.DataFrame, group_var: Union[str, List[str]]) -> pd.DataFrame:
    tx_starts = exons.groupby(group_var)['start'].min().reset_index()
    overall_start = exons['start'].min()
    
    tx_start_gaps = tx_starts.copy()
    tx_start_gaps['end'] = tx_start_gaps['start']
    tx_start_gaps['start'] = overall_start
    tx_start_gaps['seqnames'] = exons['seqnames'].iloc[0]
    tx_start_gaps['strand'] = exons['strand'].iloc[0]
    
    return tx_start_gaps

def _get_gap_map(df: pd.DataFrame, gaps: pd.DataFrame) -> dict:
    equal_hits = []
    within_hits = []

    for i, gap in gaps.iterrows():
        for j, row in df.iterrows():
            if gap['start'] == row['start'] and gap['end'] == row['end']:
                equal_hits.append({'gap_index': i, 'df_index': j})
            elif gap['start'] >= row['start'] and gap['end'] <= row['end']:
                within_hits.append({'gap_index': i, 'df_index': j})

    equal_hits = pd.DataFrame(equal_hits)
    within_hits = pd.DataFrame(within_hits)

    pure_within_hits = within_hits[~within_hits.apply(tuple, 1).isin(equal_hits.apply(tuple, 1))]

    return {'equal': equal_hits, 'pure_within': pure_within_hits}

def _get_shortened_gaps(df: pd.DataFrame, gaps: pd.DataFrame, gap_map: dict, 
                        group_var: Union[str, List[str]], target_gap_width: int) -> pd.DataFrame:
    df = df.copy()
    df['width'] = df['end'] - df['start'] + 1

    df['shorten_type'] = 'none'
    df.loc[df.index.isin(gap_map['equal']['df_index']), 'shorten_type'] = 'equal'
    df.loc[df.index.isin(gap_map['pure_within']['df_index']), 'shorten_type'] = 'pure_within'

    df['shortened_width'] = np.where(
        (df['shorten_type'] == 'equal') & (df['width'] > target_gap_width),
        target_gap_width,
        df['width']
    )

    if len(gap_map['pure_within']) > 0:
        gap_widths = gaps['end'] - gaps['start'] + 1
        sum_gap_diff = gap_map['pure_within'].groupby('df_index').apply(
            lambda x: sum(np.minimum(gap_widths.iloc[x['gap_index']], target_gap_width) - target_gap_width)
        ).reset_index(name='sum_shortened_gap_diff')

        df = df.merge(sum_gap_diff, left_index=True, right_on='df_index', how='left')
        df['shortened_width'] = np.where(
            df['shorten_type'] == 'pure_within',
            df['width'] - df['sum_shortened_gap_diff'],
            df['shortened_width']
        )
        df = df.drop(columns=['sum_shortened_gap_diff', 'df_index'])

    df = df.drop(columns=['shorten_type', 'width']).rename(columns={'shortened_width': 'width'})
    return df

def _get_rescaled_txs(exons: pd.DataFrame, introns_shortened: pd.DataFrame, 
                      tx_start_gaps_shortened: pd.DataFrame, 
                      group_var: Union[str, List[str]]) -> pd.DataFrame:
    exons = exons.copy()
    exons['width'] = exons['end'] - exons['start'] + 1

    rescaled_tx = pd.concat([exons, introns_shortened], ignore_index=True)
    rescaled_tx = rescaled_tx.sort_values(by=[group_var, 'start', 'end'] if group_var else ['start', 'end'])

    rescaled_tx['rescaled_end'] = rescaled_tx.groupby(group_var)['width'].cumsum() if group_var else rescaled_tx['width'].cumsum()
    rescaled_tx['rescaled_start'] = rescaled_tx['rescaled_end'] - rescaled_tx['width'] + 1

    if group_var is None:
        rescaled_tx['width_tx_start'] = 1
    else:
        rescaled_tx = rescaled_tx.merge(tx_start_gaps_shortened, on=group_var, suffixes=('', '_tx_start'))

    rescaled_tx['rescaled_end'] += rescaled_tx['width_tx_start']
    rescaled_tx['rescaled_start'] += rescaled_tx['width_tx_start']

    rescaled_tx['start'] = np.where(rescaled_tx['type'] == 'intron', rescaled_tx['start'] - 1, rescaled_tx['start'])
    rescaled_tx['end'] = np.where(rescaled_tx['type'] == 'intron', rescaled_tx['end'] + 1, rescaled_tx['end'])
    rescaled_tx['rescaled_start'] = np.where(rescaled_tx['type'] == 'intron', rescaled_tx['rescaled_start'] - 1, rescaled_tx['rescaled_start'])
    rescaled_tx['rescaled_end'] = np.where(rescaled_tx['type'] == 'intron', rescaled_tx['rescaled_end'] + 1, rescaled_tx['rescaled_end'])

    rescaled_tx = rescaled_tx.drop(columns=['start', 'end', 'width'])
    rescaled_tx = rescaled_tx.rename(columns={'rescaled_start': 'start', 'rescaled_end': 'end'})

    column_order = ['seqnames', 'start', 'end', 'strand'] + [col for col in rescaled_tx.columns if col not in ['seqnames', 'start', 'end', 'strand']]
    rescaled_tx = rescaled_tx[column_order]

    return rescaled_tx

def plot_transcript(exons: pd.DataFrame, introns: pd.DataFrame, 
                    group_var: str = 'transcript_name', 
                    shortened: bool = False, 
                    target_gap_width: int = 100) -> go.Figure:
    """
    Plot transcript structure using Plotly.

    Parameters:
    -----------
    exons : pd.DataFrame
        DataFrame containing exon coordinates.
    introns : pd.DataFrame
        DataFrame containing intron coordinates.
    group_var : str, optional
        Column name for grouping transcripts. Default is 'transcript_name'.
    shortened : bool, optional
        Whether to use shortened gaps. Default is False.
    target_gap_width : int, optional
        The width to shorten gaps to if shortened is True. Default is 100.

    Returns:
    --------
    go.Figure
        A Plotly figure object containing the transcript structure plot.
    """
    if shortened:
        data = shorten_gaps(exons, introns, group_var, target_gap_width)
    else:
        data = pd.concat([exons, introns], ignore_index=True)

    fig = go.Figure()

    for transcript in data[group_var].unique():
        transcript_data = data[data[group_var] == transcript]
        
        # Add exons
        exon_data = transcript_data[transcript_data['type'] == 'exon']
        fig.add_trace(go.Bar(
            x=exon_data['start'] + (exon_data['end'] - exon_data['start']) / 2,
            y=[transcript] * len(exon_data),
            width=exon_data['end'] - exon_data['start'],
            name=f"{transcript} exons",
            marker_color='blue',
            hoverinfo='text',
            hovertext=[f"Exon: {start}-{end}" for start, end in zip(exon_data['start'], exon_data['end'])]
        ))

        # Add introns
        intron_data = transcript_data[transcript_data['type'] == 'intron']
        fig.add_trace(go.Scatter(
            x=intron_data['start'],
            y=[transcript] * len(intron_data),
            mode='lines',
            line=dict(color='gray', width=2),
            name=f"{transcript} introns",
            hoverinfo='text',
            hovertext=[f"Intron: {start}-{end}" for start, end in zip(intron_data['start'], intron_data['end'])]
        ))

    fig.update_layout(
        title="Transcript Structure",
        xaxis_title="Genomic Position",
        yaxis_title="Transcript",
        barmode='overlay',
        showlegend=False
    )

    return fig

# Example usage
if __name__ == "__main__":
    # Create sample data
    exons = pd.DataFrame({
        'start': [100, 300, 600, 150, 400],
        'end': [200, 400, 700, 250, 500],
        'strand': ['+', '+', '+', '-', '-'],
        'seqnames': ['chr1'] * 5,
        'transcript_name': ['tx1', 'tx1', 'tx1', 'tx2', 'tx2'],
        'type': ['exon'] * 5
    })

    introns = pd.DataFrame({
        'start': [201, 401, 251],
        'end': [299, 599, 399],
        'strand': ['+', '+', '-'],
        'seqnames': ['chr1'] * 3,
        'transcript_name': ['tx1', 'tx1', 'tx2'],
        'type': ['intron'] * 3
    })

    # Plot original transcript structure
    #fig_original = plot_transcript(exons, introns)
    #fig_original.show()

    # Plot transcript structure with shortened gaps
    #fig_shortened = plot_transcript(exons, introns, shortened=True, target_gap_width=50)
    #fig_shortened.show()