import sys
import csv
import gzip

for file_path in sys.argv[1:]:
    chrom_pos_counts = {}
    with gzip.open(file_path, 'rt') as file:
        next(file)
        for line in file:
            if line.startswith('CHROM'):
                continue  
            fields = line.strip().split('\t')
            chrom, pos, ac_whitebritish = fields[0], fields[1], fields[4]
            if int(ac_whitebritish) > 0:
                chrom_pos = f"{chrom} {pos}"
                chrom_pos_counts[chrom_pos] = chrom_pos_counts.get(chrom_pos, 0) + 1
        
    output_file_path = f"{file_path}.counts.csv"
    with open(output_file_path, 'w', newline='') as output_file:
        writer = csv.writer(output_file)
        writer.writerow(['CHROM', 'POS', 'COUNT'])  
        for chrom_pos, count in chrom_pos_counts.items():
            chrom, pos = chrom_pos.split()
            writer.writerow([chrom, pos, count])

    print(f"Counts for file '{file_path}' written to '{output_file_path}'.")

