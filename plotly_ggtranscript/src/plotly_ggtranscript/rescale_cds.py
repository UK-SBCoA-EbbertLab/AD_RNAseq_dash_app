import pandas as pd

def rescale_cds(cds_exon_diff, gene_rescaled_exons, factor_order):
    """
    Rescales CDS regions based on exon positions and the calculated differences between CDS and exon positions.
    This function aligns the CDS regions to the rescaled exon positions and adjusts their start and end points
    accordingly.

    Parameters
    ----------
    cds_exon_diff : pd.DataFrame
        DataFrame containing the differences between the start and end positions of exons and CDS regions.
        Expected columns:
        - 'd_start': Absolute difference between exon start and CDS start positions.
        - 'd_end': Absolute difference between exon end and CDS end positions.
        - Any columns necessary for joining (e.g., 'transcript_id', 'gene_id').
    gene_rescaled_exons : pd.DataFrame
        DataFrame containing the rescaled exon positions.
        Expected columns:
        - 'start': Rescaled start position of the exon.
        - 'end': Rescaled end position of the exon.
        - Any columns necessary for joining.
    factor_order : list
        List specifying the desired order of 'transcript_id' for categorical ordering.

    Returns
    -------
    gene_rescaled_cds : pd.DataFrame
        DataFrame containing the rescaled CDS positions, with adjusted 'start' and 'end' positions,
        and 'transcript_id' ordered according to 'factor_order'.
    """

    # Step 1: Prepare the CDS DataFrame
    # - Assign a new column 'type' with the value "CDS"
    # - Drop unnecessary columns: 'c_start', 'c_end', 'e_start', 'e_end' if they exist

    cds_prepared = (
        cds_exon_diff
        .assign(type="CDS")
        .drop(columns=['c_start', 'c_end', 'e_start', 'e_end'], errors='ignore')
    )

    # Debug statement to check the prepared CDS DataFrame
    print("Prepared CDS DataFrame:\n", cds_prepared.head())

    # Step 2: Prepare the Exon DataFrame
    # - Rename 'start' to 'e_start' and 'end' to 'e_end' to avoid column name conflicts
    # - Drop the 'type' column if it exists

    exons_prepared = (
        gene_rescaled_exons
        .rename(columns={'start': 'e_start', 'end': 'e_end'})
        .drop(columns=['type'], errors='ignore')
    )

    # Debug statement to check the prepared Exon DataFrame
    print("Prepared Exon DataFrame:\n", exons_prepared.head())

    # Step 3: Identify common columns for joining
    # - Find columns that are present in both DataFrames to use as keys for the join

    common_columns = [col for col in cds_prepared.columns if col in exons_prepared.columns]
    if not common_columns:
        raise ValueError("No common columns to perform join on. Ensure both DataFrames have common keys.")

    print("Common columns for joining:", common_columns)

    # Step 4: Perform the left join on the common columns
    # - This aligns the CDS data with the corresponding rescaled exon positions

    gene_rescaled_cds = pd.merge(
        cds_prepared,
        exons_prepared,
        how='left',
        on=common_columns
    )

    # Debug statement to check the merged DataFrame
    print("DataFrame after left join:\n", gene_rescaled_cds.head())

    # Step 5: Calculate the adjusted 'start' and 'end' positions for the CDS regions
    # - 'start' is adjusted by adding 'd_start' to 'e_start'
    # - 'end' is adjusted by subtracting 'd_end' from 'e_end'

    gene_rescaled_cds['start'] = gene_rescaled_cds['e_start'] + gene_rescaled_cds['d_start']
    gene_rescaled_cds['end'] = gene_rescaled_cds['e_end'] - gene_rescaled_cds['d_end']

    # Debug statement to check the DataFrame after adjusting positions
    print("DataFrame after adjusting 'start' and 'end':\n", gene_rescaled_cds[['start', 'end']].head())

    # Step 6: Drop unnecessary columns used for calculations
    # - Remove 'e_start', 'e_end', 'd_start', 'd_end' as they are no longer needed

    gene_rescaled_cds = gene_rescaled_cds.drop(
        columns=['e_start', 'e_end', 'd_start', 'd_end'],
        errors='ignore'
    )

    # Debug statement to check the final DataFrame before setting 'transcript_id' as categorical
    print("Final DataFrame before setting 'transcript_id' as categorical:\n", gene_rescaled_cds.head())

    # Step 7: Set 'transcript_id' as a categorical variable with the specified order
    # - This ensures that the 'transcript_id' column follows the desired ordering

    if 'transcript_id' in gene_rescaled_cds.columns:
        gene_rescaled_cds['transcript_id'] = pd.Categorical(
            gene_rescaled_cds['transcript_id'],
            categories=factor_order,
            ordered=True
        )
    else:
        raise KeyError("'transcript_id' column not found in the DataFrame.")

    # Debug statement to confirm 'transcript_id' is set as categorical
    print("DataFrame after setting 'transcript_id' as categorical:\n", gene_rescaled_cds.head())

    return gene_rescaled_cds
