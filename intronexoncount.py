 """
    Parses a GTF file and extracts more intron and exon stuff
    spits out inton exon counts
"""

import os
import re
import pandas as pd
from collections import defaultdict

def parse_gtf(file_path):
    """Parses a GTF file and returns a dictionary of gene_id to exon counts."""
    gene_exons = defaultdict(int)
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith("#"):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != "exon":
                    continue

                attributes = fields[8]
                gene_id_match = re.search(r'gene_id "(.*?)"', attributes)
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    gene_exons[gene_id] += 1
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
    
    return gene_exons

def calculate_averages(gene_exons):
    """Calculates the average number of exons and introns per gene."""
    total_exons = sum(gene_exons.values())
    total_genes = len(gene_exons)

    avg_exons = total_exons / total_genes if total_genes > 0 else 0
    avg_introns = (total_exons - total_genes) / total_genes if total_genes > 0 else 0

    return avg_exons, avg_introns

def main():
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

    results = []

    for gtf_file in gtf_files:
        if not os.path.exists(gtf_file):
            print(f"File not found: {gtf_file}")
            continue

        gene_exons = parse_gtf(gtf_file)
        avg_exons, avg_introns = calculate_averages(gene_exons)

        results.append({
            "File": gtf_file,
            "Average Exons per Gene": avg_exons,
            "Average Introns per Gene": avg_introns
        })

    # Convert results to a DataFrame for better visualization
    results_df = pd.DataFrame(results)
    print(results_df)

if __name__ == "__main__":
    main()
