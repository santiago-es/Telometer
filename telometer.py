#!/usr/bin/env python3
# Telometer v0.51
# Created by: Santiago E Sanchez
# Artandi Lab, Stanford University, 2023
# Measures telomeres from ONT or PacBio long reads aligned to a T2T genome assembly
# Simple Usage: python telometer.py -b sorted_t2t.bam -o output.tsv

import pysam
import re
import csv
import argparse
import pandas as pd

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement[base] for base in reversed(seq))

def calculate_telomere_length():
    # required inputs: bam_file_path, output_file_path, chemistry
    parser = argparse.ArgumentParser(description='Calculate telomere length from a BAM file.')
    parser.add_argument('-b', '--bam', help='The path to the sorted BAM file.', required=True)
    parser.add_argument('-o', '--output', help='The path to the output file.', required=True)
    parser.add_argument('-c', '--chemistry', default="r9", help="Sequencing chemistry (r9 or r10). Optional", required=False)
    parser.add_argument('-m', '--minreadlen', default=1000, help='Minimum read length to consider (Default: 1000). Optional', required=False)
    args = parser.parse_args()
    bam_file = pysam.AlignmentFile(args.bam, "rb")


    if args.chemistry == 'r10':
        adapters = ['TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT', 'GCAATACGTAACTGAACGAAGTACAGG']
        adapters_rc = [reverse_complement(adapt) for adapt in adapters]
    else:
        adapters = ['TTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT', 'GCAATACGTAACTGAACGAAGT']
        adapters_rc = [reverse_complement(adapt) for adapt in adapters]
    for y in adapters_rc:
        adapters.append(y)
    telomere_repeats = ['GGCCA','CCCTAA', 'TTAGGG', 'CCCTGG', 'CTTCTT', 'TTAAAA', 'CCTGG']
    telomere_repeats_rc = [reverse_complement(repeat) for repeat in telomere_repeats]
    for x in telomere_repeats_rc:
        telomere_repeats.append(x)

    telomere_repeats_re = "|".join(f'({repeat}){{2,}}' for repeat in telomere_repeats)

    # Extract length of the reference sequence from the BAM file header, for debugging
    # reference_genome_length = bam_file.header['SQ'][0]['LN']


    highest_mapping_quality = {}


    results = []
    p_count = 0
    q_count = 0
    rev_count = 0
    fwd_count = 0
    p_tel = 0
    q_tel = 0
    for read in bam_file:
        alignment_start = read.reference_start
        alignment_end = read.reference_end
        seq = read.query_sequence
        if read.is_unmapped:
            continue
        reference_genome_length = bam_file.get_reference_length(read.reference_name)
        if "M" in read.reference_name or len(seq) < args.minreadlen:
            continue
        if read.is_reverse:
            rev_count += 1
            direction = "rev"
            seq = reverse_complement(seq)

        else:
            fwd_count += 1
            direction = "fwd"
        if (alignment_start < 15000 or alignment_start > reference_genome_length - 30000):
            # find start of telomeric repeat region
            if alignment_start < 15000 and "q" not in str(read.reference_name):
                arm = "p"
                p_count += 1
            elif alignment_start < 15000 and "q" in str(read.reference_name):
                arm = "q"
                q_count += 1
            else:
                arm = "q"
                q_count += 1

            telomere_start = [m.start() for m in re.finditer(telomere_repeats_re, seq)]
            if telomere_start:
                telomere_start = telomere_start[0]
                ## Terminal telomere check
                if telomere_start > 100 and (len(seq) - telomere_start > 200):
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
                    if arm == "p":
                        p_tel += 1
                    else:
                        q_tel += 1
                    highest_mapping_quality[read.query_name] = read.mapping_quality
                    results.append({
                        'chromosome': read.reference_name,
                        'reference_start': alignment_start,
                        'reference_end': alignment_end,
                        'telomere_length': telomere_length,
                        'read_id': read.query_name,
                        'mapping_quality': read.mapping_quality,
                        'read_length': len(seq),
                        'arm': arm,
                        'direction': direction
                    })

    bam_file.close()

    # Write the results to the output file
    with open(args.output, 'w') as output_file:
    # with open(output_file_path, 'w', newline='') as output_file:
        writer = csv.DictWriter(output_file, fieldnames=results[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(results)

    #print("p arm seen "+str(p_count))
    #print("q arm seen "+str(q_count))
    #print("fwd "+str(fwd_count))
    #print("rev "+str(rev_count))
    #print("p tel "+ str(p_tel))
    #print("q tel "+str(q_tel))

if __name__ == "__main__":
    calculate_telomere_length()
