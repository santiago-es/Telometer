# telomap-ont
A simple regular expression based method for measuring telomere length from noisy Oxford Nanopore Technology reads


Dependencies: pysam, regex, csv, argparse

Simple Usage: python3 telomap-ont.py -b /path/to/sorted.bam -o /path/to/output.tsv

FASTQ reads should be aligned to the T2T-CHM13-2.0 genome with minimap2 or winnowmap.

The logic of this script is based on telomap for PacBio (see: https://github.com/cytham/telomap/) and [CY Tham et al Nat Comms 2023](https://www-nature-com.laneproxy.stanford.edu/articles/s41467-023-35823-7#Sec1)

Namely, it looks for telomeric repeats using regular expression itertools and measures telomeres from the sequencing adapter sequence to the last telomeric repeat of the form 5'-TTAGGG-3' or 5'-AATCCC-3'.
Because ONT reads are noisy and frequently miscall telomeres (see [Tan et al Gen. Bio. 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02751-6)) in stereotypical modes, my script also counts these frequently miscalled motifs as telomeric repeats.

Additionally, my script only searches reads which align to the first or last 15,000 bp of their reference chromosome and only measures telomeres from reads longer than 5kb to ensure any analyzed read would be sufficiently long to contain likely intact telomres. It then checks that the first or last 100 bp of a read are telomere-rich to ensure telomere measurements are from terminal and not interstitial telomere sequences.
