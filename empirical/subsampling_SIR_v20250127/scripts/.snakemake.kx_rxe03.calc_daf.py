
######## Snakemake header ########
import sys; sys.path.insert(0, "/Users/margaretsteiner/miniconda3/envs/snakemake/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X \x00\x00\x00metadata/ancestral_vars_chr1.txtq\x06XK\x00\x00\x00results/freq/chr1_missense_centerE9N9geo150000_nSIR10000_nSIRreps10.SIRfreqq\x07e}q\x08(X\x06\x00\x00\x00_namesq\t}q\n(X\x16\x00\x00\x00anc_alleles_file_cleanq\x0bK\x00N\x86q\x0cX\x0c\x00\x00\x00frq_file_snvq\rK\x01N\x86q\x0euh\x0bh\x06h\rh\x07ubX\x06\x00\x00\x00outputq\x0fcsnakemake.io\nOutputFiles\nq\x10)\x81q\x11XF\x00\x00\x00results/daf/chr1_missense_centerE9N9geo150000_nSIR10000_nSIRreps10.dafq\x12a}q\x13(h\t}q\x14X\x07\x00\x00\x00daf_snvq\x15K\x00N\x86q\x16sh\x15h\x12ubX\x06\x00\x00\x00paramsq\x17csnakemake.io\nParams\nq\x18)\x81q\x19}q\x1ah\t}q\x1bsbX\t\x00\x00\x00wildcardsq\x1ccsnakemake.io\nWildcards\nq\x1d)\x81q\x1e(X\x01\x00\x00\x001q\x1fX\x08\x00\x00\x00missenseq X\x13\x00\x00\x00centerE9N9geo150000q!X\x05\x00\x00\x0010000q"X\x02\x00\x00\x0010q#e}q$(h\t}q%(X\x05\x00\x00\x00chromq&K\x00N\x86q\'X\x08\x00\x00\x00varclassq(K\x01N\x86q)X\x08\x00\x00\x00scenarioq*K\x02N\x86q+X\x05\x00\x00\x00nsampq,K\x03N\x86q-X\x05\x00\x00\x00nrepsq.K\x04N\x86q/uX\x05\x00\x00\x00chromq0h\x1fX\x08\x00\x00\x00varclassq1h X\x08\x00\x00\x00scenarioq2h!X\x05\x00\x00\x00nsampq3h"X\x05\x00\x00\x00nrepsq4h#ubX\x07\x00\x00\x00threadsq5K\x01X\t\x00\x00\x00resourcesq6csnakemake.io\nResources\nq7)\x81q8(K\x01K\x01e}q9(h\t}q:(X\x06\x00\x00\x00_coresq;K\x00N\x86q<X\x06\x00\x00\x00_nodesq=K\x01N\x86q>uh;K\x01h=K\x01ubX\x03\x00\x00\x00logq?csnakemake.io\nLog\nq@)\x81qA}qBh\t}qCsbX\x06\x00\x00\x00configqD}qEX\x04\x00\x00\x00ruleqFX\x13\x00\x00\x00derived_allele_freqqGub.')
######## Original script #########
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
