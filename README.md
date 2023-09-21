# Telometer

v0.5
A simple regular expression based method for measuring telomere length from long read sequencing

Dependencies: pysam, regex, pandas

Simple Usage: 
```
python3 telometer.py -b /path/to/sorted.bam -o /path/to/output.tsv
```

# Workflow

FASTQ reads should be aligned to the T2T-CHM13-2.0 genome with minimap2 or winnowmap. This repository contains a useful reference (t2t-and-subtel.fa) which combines the hs1 (CHM13 v2.0) assembly with alternative subtelomere assemblies ([Stong 2014](https://pubmed.ncbi.nlm.nih.gov/24676094/)). 

```
minimap2 -ax map-ont \
 -t 16 \ 
-N 5 \ -Y 
-L \ 
/path/to/reference/t2t-and-subtel.fa \
/path/to/fastq_dir/*.fastq  \
-o output.sam
```
Convert the sam output to bam, sort, and index

```
samtools view -bho output.bam output.sam
samtools sort -o output-sort.bam output.bam
samtools index output-sort.bam
```
Run Telometer

```
python3 telometer.py -b /path/to/output-sort.bam -o /path/to/output.tsv
```


Telometer looks for telomeric repeats using regular expressions and measures telomeres from the sequencing adapter sequence to the last telomeric repeat of the form 5'-TTAGGG-3' or 5'-AATCCC-3'.
Because ONT reads are noisy and frequently miscall telomeres (see [Tan et al Gen. Bio. 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02751-6)) in stereotypical modes, Telometer also counts these frequently miscalled motifs as telomeric repeats.

Additionally, my script only searches reads which align to the first or last several thousand base pairs of their reference chromosome and only measures telomeres from reads longer than 1000 bp to ensure any analyzed read would be sufficiently long to contain likely intact telomres. It then checks that the first or last 100 bp of a read are telomere-rich to ensure telomere measurements are from terminal and not interstitial telomere sequences.

If this script is helpful please cite the original article, once it is published. 
