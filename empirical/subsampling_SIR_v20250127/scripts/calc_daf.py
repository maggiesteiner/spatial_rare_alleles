import sys


def load_ancestral_alleles(anc_file):
    """Load ancestral alleles into a dictionary for fast lookup."""
    anc_dict = {}
    count = 0
    with open(anc_file, 'r') as f:
        # Skip the header line
        f.readline()
        for line in f:
            chrom, pos, allele = line.strip().split("\t")  # Split on tab for the ancestral file
            anc_dict[(chrom, pos)] = allele.upper()  # Store in uppercase
            count += 1
            if count % 500000 == 0:  # Print every 500K records
                print(f"Loaded {count} ancestral alleles...", file=sys.stderr)
    print(f"Finished loading {count} ancestral alleles.", file=sys.stderr)
    return anc_dict  # O(1) lookup time


def process_snps(anc_file, snp_file, output_file):
    """Process SNP file and compute derived allele frequencies."""

    # Load ancestral alleles
    ancestral = load_ancestral_alleles(anc_file)

    processed_count = 0
    skipped_count = 0

    with open(snp_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write header
        outfile.write("ITER\tCHR\tPOS\tA1\tA2\tMAF\tDerived_Frequency\tNOBS\n")

        # Skip the header line of SNP file (assuming the first line is a header)
        infile.readline()

        for line in infile:
            parts = line.strip().split()  # Split on any whitespace (spaces or tabs)

            # Skip lines that don't have exactly 7 columns
            if len(parts) != 7:
                print(f"Skipping malformed line: {line.strip()}", file=sys.stderr)
                skipped_count += 1
                continue  # Skip malformed line

            iter, chr, pos, maf, nobs, a1, a2 = parts[:7]

            # Convert to uppercase for case-insensitive comparison
            a1, a2 = a1.upper(), a2.upper()
            ancestral_allele = ancestral.get((chr, pos), None)

            # Skip missing or invalid alleles
            if ancestral_allele in {None, ".", "-"} or a1 == "." or a2 == ".":
                print(f"Skipping SNP due to missing/invalid allele: {chr} {pos}", file=sys.stderr)
                skipped_count += 1
                continue  # Skip this SNP

            # Determine derived frequency
            if a1 == ancestral_allele:
                derived_freq = 1 - float(maf)
            elif a2 == ancestral_allele:
                derived_freq = float(maf)
            else:
                print(
                    f"Neither A1 nor A2 matches ancestral for: {chr} {pos} (Ancestral: {ancestral_allele}, A1: {a1}, A2: {a2})",
                    file=sys.stderr)
                skipped_count += 1
                continue  # Skip SNP if neither allele matches ancestral

            # Write output
            outfile.write(f"{iter}\t{chr}\t{pos}\t{a1}\t{a2}\t{maf}\t{derived_freq}\t{nobs}\n")
            processed_count += 1

            # Print progress updates every 1M processed SNPs
            if processed_count % 1000000 == 0:
                print(f"Processed {processed_count} SNPs...", file=sys.stderr)

    print(f"Finished processing SNPs: {processed_count} written, {skipped_count} skipped.", file=sys.stderr)


# Snakemake inputs/outputs
process_snps(snakemake.input.anc_alleles_file_clean, snakemake.input.frq_file_snv, snakemake.output.daf_snv)
