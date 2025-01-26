"""
    Parses a GTF file and extracts more intron and exon stuff
    spits out inton exon counts
"""
import pandas as pd
import statistics

def parse_gtf(gtf_file):
    gene_lengths = []

    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue  # Skip header lines
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # Skip lines that do not have enough fields (e.g., empty lines)
            
            # Check if the feature type is "gene"
            if "gene" in fields[2]:
                start = int(fields[3])  # Gene start position
                end = int(fields[4])    # Gene end position
                length = end - start + 1  # Calculate gene length

                gene_lengths.append(length)  # Append gene length to list

    return gene_lengths

def calculate_statistics(lengths):

    mean_length = sum(lengths) / len(lengths) if lengths else 0
    median_length = statistics.median(lengths) if lengths else 0
    return mean_length, median_length

def main():
    # List of GTF file names
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
    
    # Process each file
    for gtf_file in gtf_files:
        print(f"Processing file: {gtf_file}")
        try:
            gene_lengths = parse_gtf(gtf_file)
            if not gene_lengths:
                print(f"No gene features found in {gtf_file}.")
                continue

            mean_length, median_length = calculate_statistics(gene_lengths)
            print(f"File: {gtf_file}")
            print(f"Mean gene length: {mean_length}")
            print(f"Median gene length: {median_length}\n")
        except FileNotFoundError:
            print(f"File not found: {gtf_file}")
        except Exception as e:
            print(f"An error occurred while processing {gtf_file}: {e}")

if __name__ == "__main__":
    main()
