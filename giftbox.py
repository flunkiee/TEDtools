import os
import re
import pandas as pd
from collections import defaultdict
from tqdm import tqdm

def parse_gtf(file_path):
    """Parse a GTF file and return a DataFrame and a dictionary of gene_id to exon counts."""
    gene_exons = defaultdict(int)
    data = []

    try:
        with open(file_path, 'r') as file:
            for line in tqdm(file, desc=f"Processing {os.path.basename(file_path)}", unit=" lines"):
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                feature = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]
                gene_id_match = re.search(r'gene_id "(.*?)"', attributes)
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    if feature == "exon":
                        gene_exons[gene_id] += 1
                    data.append((fields[0], fields[1], feature, start, end, fields[5], fields[6], fields[7], gene_id))
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return pd.DataFrame(), {}

    df = pd.DataFrame(data, columns=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id'])
    return df, gene_exons

def calculate_averages(gene_exons):
    """Calculate the average number of exons and introns per gene."""
    total_exons = sum(gene_exons.values())
    total_genes = len(gene_exons)
    avg_exons = total_exons / total_genes if total_genes > 0 else 0
    avg_introns = (total_exons - total_genes) / total_genes if total_genes > 0 else 0
    return avg_exons, avg_introns

def calculate_exon_intron_lengths(df):
    """Calculate mean and median exon and intron lengths."""
    exons = df[df['feature'] == 'exon'].copy()
    if exons.empty:
        return 0, 0, 0, 0
    exons['length'] = exons['end'] - exons['start'] + 1
    exons_sorted = exons.sort_values(['gene_id', 'start'])
    introns = []
    for gene_id, group in exons_sorted.groupby('gene_id'):
        for i in range(1, len(group)):
            intron_length = group.iloc[i]['start'] - group.iloc[i-1]['end'] - 1
            if intron_length > 0:
                introns.append(intron_length)
    mean_exon_length = exons_sorted['length'].mean()
    median_exon_length = exons_sorted['length'].median()
    mean_intron_length = pd.Series(introns).mean() if introns else 0
    median_intron_length = pd.Series(introns).median() if introns else 0
    return mean_exon_length, median_exon_length, mean_intron_length, median_intron_length

def main():
    # Prompt user to input GTF files
    gtf_files = input("Enter the GTF files to process, separated by commas: ").split(',')

    results = []

    for gtf_file in gtf_files:
        gtf_file = gtf_file.strip()
        if not os.path.exists(gtf_file):
            print(f"File not found: {gtf_file}")
            continue

        df, gene_exons = parse_gtf(gtf_file)
        if df.empty:
            print(f"No valid data found in file: {gtf_file}")
            continue

        avg_exons, avg_introns = calculate_averages(gene_exons)
        mean_exon_length, median_exon_length, mean_intron_length, median_intron_length = calculate_exon_intron_lengths(df)

        results.append({
            "File": gtf_file,
            "Average Exons per Gene": avg_exons,
            "Average Introns per Gene": avg_introns,
            "Mean Exon Length": mean_exon_length,
            "Median Exon Length": median_exon_length,
            "Mean Intron Length": mean_intron_length,
            "Median Intron Length": median_intron_length
        })

    results_df = pd.DataFrame(results)
    print(results_df)

if __name__ == "__main__":
    main()