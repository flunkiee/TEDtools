 """
    parses GTF file, intsona and exons

    
    spits out : A list of intron and exon stuff
    """

import pandas as pd

# List of GTF files
gtf_files = [
    "A.australiensis_braker_EDTA_B3.gtf",
    "C_formosanus_EDTAbraker.gtf",
    "C_fukii_EDTAbraker.gtf",
    "G_aquaticus_EDTAbraker.gtf",
    "G_montsenyensis_EDTAbraker(1).gtf",
    "G_parenensis_EDTAbraker(1).gtf",
    "N.muniae_EDTAbraker3.gtf",
    "pvarius braker.gtf"
]

# Open a file to write the output
with open('gtf_analysis_output.txt', 'w') as output_file:
    # Iterate over each GTF file and calculate lengths
    for gtf_file in gtf_files:
        try:
            # Load the GTF file into a pandas DataFrame
            df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

            # Set column names
            df.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

            # Filter exons and create a new DataFrame for further processing
            exons = df[df['feature'] == 'exon'].copy()

            # Calculate exon lengths and assign them directly to the 'length' column
            exons['length'] = exons['end'] - exons['start'] + 1

            # Extract gene_id and assign it to a new 'gene_id' column
            exons['gene_id'] = exons['attribute'].str.extract('gene_id "([^"]+)"')

            # Sort exons by start position within each gene
            exons_sorted = exons.sort_values(['gene_id', 'start'])

            # Calculate intron lengths
            introns = []
            for gene_id, group in exons_sorted.groupby('gene_id'):
                for i in range(1, len(group)):
                    intron_length = group.iloc[i]['start'] - group.iloc[i-1]['end'] - 1
                    if intron_length > 0:
                        introns.append(intron_length)

            # Calculate mean and median lengths
            mean_exon_length = exons_sorted['length'].mean()
            median_exon_length = exons_sorted['length'].median()
            mean_intron_length = pd.Series(introns).mean()
            median_intron_length = pd.Series(introns).median()

            # Write the results to the file
            output_file.write(f"Results for {gtf_file}:\n")
            output_file.write(f"  Mean exon length: {mean_exon_length}\n")
            output_file.write(f"  Median exon length: {median_exon_length}\n")
            output_file.write(f"  Mean intron length: {mean_intron_length}\n")
            output_file.write(f"  Median intron length: {median_intron_length}\n")
            output_file.write("\n" + "-"*50 + "\n")

        except Exception as e:
            output_file.write(f"Error processing {gtf_file}: {e}\n")
            output_file.write("\n" + "-"*50 + "\n")

