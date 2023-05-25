import pysam
import re
import csv



def calculate_telomere_length(bam_file_path, output_file_path):
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")

    adapters = ['TTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT', 'GCAATACGTAACTGAACGAAGT']
    telomere_repeats = ['CCCTAA', 'TTAGGG', 'CCCTGG', 'CTTCTT', 'TTAAAA', 'CCTGG']

    telomere_repeats_re = "|".join(f'({repeat}){{2,}}' for repeat in telomere_repeats)

    # Extract length of the reference sequence from the BAM file header
    reference_genome_length = bam_file.header['SQ'][0]['LN']

    # Holds the highest mapping quality found for each read
    highest_mapping_quality = {}

    # The results to be written to the output file
    results = []

    for read in bam_file:
        alignment_start = read.reference_start
        alignment_end = read.reference_end
        seq = read.query_sequence

        # Only consider reads with length greater than 3000
        if len(seq) > 3000 and (alignment_start < 15000 or alignment_end > reference_genome_length - 15000):
            # find start of telomeric repeat region
            telomere_start = [m.start() for m in re.finditer(telomere_repeats_re, seq)]

            if telomere_start:
                telomere_start = telomere_start[0]

                # Terminal telomere check
                if telomere_start > 100 and (len(seq) - telomere_start > 100):
                    continue

                # find the end of the telomeric repeat region by locating the adapter sequence
                telomere_end = min((seq.find(adapter) for adapter in adapters), default=-1)
                if telomere_end == -1:
                    telomere_end = len(seq)

                telomere_region = seq[telomere_start:telomere_end]

                # get all telomeric repeat occurrences in the telomere_region
                telomere_repeat = [m.group() for m in re.finditer('|'.join(telomere_repeats), telomere_region)]

                telomere_length = len(''.join(telomere_repeat))

                # Only record this measurement if it is the highest quality for this read ID
                if read.query_name not in highest_mapping_quality or read.mapping_quality > highest_mapping_quality[read.query_name]:
                    highest_mapping_quality[read.query_name] = read.mapping_quality
                    results.append({
                        'chromosome': read.reference_name,
                        'reference_start': alignment_start,
                        'reference_end': alignment_end,
                        'telomere_length': telomere_length,
                        'read_id': read.query_name,
                        'mapping_quality': read.mapping_quality,
                        'read_length': len(seq)
                    })

    bam_file.close()

    # Write the results to the output file
    with open(output_file_path, 'w', newline='') as output_file:
        writer = csv.DictWriter(output_file, fieldnames=results[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(results)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate telomere length from a BAM file.')
    parser.add_argument('-b', '--bam', help='The path to the sorted BAM file.', required=True)
    parser.add_argument('-o', '--output', help='The path to the output file.', required=True)

    args = parser.parse_args()

    calculate_telomere_length(args.bam, args.output)
