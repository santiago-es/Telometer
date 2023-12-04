# Telometer

![Telometer Logo](https://i.imgur.com/te0QfrR.png)

v0.5
A simple regular expression based method for measuring telomere length from long read sequencing

Dependencies: pysam, pandas

Simple Usage: 
```
telometer -b /path/to/sorted.bam -o /path/to/output.tsv
```

# Workflow

Install telometer
```
conda install -c bioconda telometer
```

Download the hs1 assembly from https://github.com/marbl/CHM13 (chm13v2.0.fa)

Append [Stong 2014](https://pubmed.ncbi.nlm.nih.gov/24676094/) subtelomere assemblies and index:
```
cat chm13v2.0.fa stong_subtels.fa > t2t-and-subtel.fa
samtools faidx t2t-and-subtel.fa
```

FASTQ reads should be aligned to the T2T-CHM13-2.0 genome with minimap2 or winnowmap.   

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

To save space, it is recommended to delete the initial sam and unsorted bam outputs and compress sorted bams. 

Run Telometer

```
telometer -b /path/to/output-sort.bam -o /path/to/output.tsv -m [minimum read length to consider, default 1000] -c [chemistry (r9 or r10), default r9]
```

Telometer looks for telomeric repeats using regular expressions and measures telomeres from the sequencing adapter sequence to the last telomeric repeat of the form 5'-TTAGGG-3' or 5'-AATCCC-3'.
Because ONT reads are noisy and frequently miscall telomeres (see [Tan et al Gen. Bio. 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02751-6)) in stereotypical modes, Telometer also counts these frequently miscalled motifs as telomeric repeats. That said, since this code was created improved telomere basecalling has been integrated into the default R10 chemistry dorado basecalling model and R10 high accuracy basecalling with dorado is now the recommended sequencing chemistry and basecalling model for telomere measurement. 

Additionally, my script only searches reads which align to the first or last several thousand base pairs of their reference chromosome and only measures telomeres from reads longer than 1000 bp to ensure any analyzed read would be sufficiently long to contain likely intact telomres. It then checks that the first or last 100 bp of a read are telomere-rich to ensure telomere measurements are from terminal and not interstitial telomere sequences.

By default, telometer only considers reads with read length greater than 1000 bp and this minimum is recommended for telomere capture libraries. For whole genome sequencing, this should be raised to at least 3000 bp.

A benchtop protocol for performing telomere captre library preparation in simplex or multiplex, please see TelometerLibraryPrep.docx in this repo. 

If this script is helpful please cite the original article:

Sanchez, S. E. et al. Digital telomere measurement by long-read sequencing distinguishes healthy aging from disease. _bioRxiv_ 2023. 

[Preprint at https://doi.org/10.1101/2023.11.29.569263]((https://www.biorxiv.org/content/10.1101/2023.11.29.569263v1)https://www.biorxiv.org/content/10.1101/2023.11.29.569263v1)


