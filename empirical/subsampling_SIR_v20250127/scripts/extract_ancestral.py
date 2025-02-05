# extract_ancestral.py

# Open the FASTA file and output file using Snakemake I/O
def extract_ancestral_alleles(fasta_file, output_file,chrom=1):
    # Initialize variables
    sequences = {}
    current_chromosome = ""
    sequence = ""

    # Open the FASTA file and read it line by line
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            line = line.strip()

            # If the line starts with '>', it's a header line (chromosome info)
            if line.startswith('>'):
                # If we have a previous sequence, store it
                if current_chromosome:
                    sequences[current_chromosome] = sequence
                # Extract the chromosome name (e.g., 'chr1')
                current_chromosome = chrom#line.split(":")[2]  # Extract chromosome part from the header (e.g., 'chr1')
                sequence = ""  # Reset sequence for the next chromosome
            else:
                # Otherwise, this is part of the sequence for the current chromosome
                sequence += line

        # Don't forget to store the last chromosome's sequence
        if current_chromosome:
            sequences[current_chromosome] = sequence

    # Open the output file to write the results
    with open(output_file, 'w') as output:
        output.write("CHR\tPOS\tANCESTRAL\n")

        # Process the entire sequence for each chromosome
        for chr_name, seq in sequences.items():
            for pos in range(len(seq)):
                # The position in the FASTA file is 1-based, so add 1 to pos
                ancestral_allele = seq[pos]  # Assume the ancestral allele is the base at this position

                # Write the result to the output file
                output.write(f"{chr_name}\t{pos + 1}\t{ancestral_allele}\n")
                print(f"{chr_name}\t{pos + 1}\t{ancestral_allele}\n")
    print(f"Ancestral output has been written to: {output_file}")

# Call the function with Snakemake I/O
if __name__ == "__main__":
    # import snakemake

    fasta_file = snakemake.input[0]  # Get the input FASTA file from Snakemake
    output_file = snakemake.output[0]  # Get the output file from Snakemake

    extract_ancestral_alleles(fasta_file, output_file)
