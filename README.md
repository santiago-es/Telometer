# Telometer

![Telometer Logo](https://i.imgur.com/te0QfrR.png)

v1.1
A simple regular expression based method for measuring telomere length from long read sequencing

Dependencies: pysam, pandas, regex, samtools, minimap2, scipy (and their associated dependencies)

Simple Usage: 
```
pip install telometer==1.1
telometer -b /path/to/sorted.bam -o /path/to/output.tsv
```
# Description

Telometer looks for telomeric repeats using regular expressions and measures telomeres from the sequencing adapter sequence to the last telomeric repeat of the form 5'-TTAGGG-3' or 5'-AATCCC-3'.
Because ONT reads are noisy and frequently miscall telomeres (see [Tan et al Gen. Bio. 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02751-6)) in stereotypical modes, Telometer originally counted these frequently miscalled motifs as telomeric repeats. That said, since this code was initially created improved telomere basecalling has been integrated into the default R10 chemistry dorado basecalling model and R10 high accuracy basecalling with dorado is now the recommended sequencing chemistry and basecalling model for telomere measurement and support for R9 has been deprecated.

Additionally, Telometer only searches reads which align to the first or last several thousand base pairs of their reference chromosome and only measures telomeres from reads longer than 1000 bp to ensure any analyzed read would be sufficiently long to contain likely intact telomeres. It then checks that the first or last 100 bp of a read are telomere-rich to ensure telomere measurements are from terminal and not interstitial telomere sequences.

By default, telometer only considers reads with read length greater than 1000 bp and this minimum is recommended for telomere capture libraries. For whole genome sequencing, this should be raised to 4000 bp.

The sub/telomeric boundary in humans tends to contain stretches of highly variable length which consist of canonical telomere motifs with one mismatch (see [Stephens and Kocher, 2024](https://link.springer.com/article/10.1186/s12859-024-05807-5)). Occasionally, these stretches are internal to a longer stretch of canonical telomere motifs. The latest version of telometer accounts for internal stretches of variants with 1 bp mismatch and counts them within the telomere length, but this may change as we learn more about telomere measurement by long read sequencing. 

For a benchtop protocol for performing telomere capture library preparation in simplex or multiplex, please see TelometerLibraryPrep.docx in this repo. 

If this code or library prep method is helpful, please cite the original article:

[Sanchez, S. E. et al. Digital telomere measurement by long-read sequencing distinguishes healthy aging from disease. _Nature Communications_ 2024](https://www.nature.com/articles/s41467-024-49007-4)



# Output Structure

| Column   | Description |
| -------- | ------- |
| chromosome | Identity of Aligned Contig for Read (As listed in reference)    |
| reference_start / *_end | Start and End of Alignment in Reference     |
| telomere_length    | Telomere Length    |
| read_id| read_id assigned by sequencer |
| read_length | read length |
| mapping_quality| MAPQ score given by minimap2 |
| arm | p or q chromosome arm |
| direction| forward (fwd) or reverse (rev) read orientation |

# Workflow

Install telometer
```
pip install telometer==0.95
```

Download the latest human t2t assembly from https://github.com/marbl/CHM13 (chm13v2.0.fa)

Append [Stong 2014](https://pubmed.ncbi.nlm.nih.gov/24676094/) subtelomere assemblies and index (stong_subtels.fa provided with source):
```
cat chm13v2.0.fa stong_subtels.fa > t2t-and-subtel.fa
samtools faidx t2t-and-subtel.fa
```

FASTQ reads should be aligned to the Stong+T2T-CHM13-2.0 genome with minimap2.   

```
minimap2 -ax map-ont \
 -t [Max # of Threads for Your Machine] \ 
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
telometer -b /path/to/output-sort.bam -o /path/to/output.tsv -m [minimum read length to consider, default 1000 for telomere capture experiments, use 4000 for WGS] -g [maximum non-telomeric gap tolerated within telomere sequence, default 100]
```
Minimal test data subsampled from a telomere capture experiment is included with source. To test:

```
telometer -b minimal_tels.bam -o output.tsv 
```
The minimal dataset should produce measurements with the following summary statistics using telometer default settings: 
```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    583    2714    3767    3861    4711   20233
```

Other options: 

| Argument      | Short Flag | Type    | Description                                                                                  |
| ------------- | ---------- | ------- | -------------------------------------------------------------------------------------------- |
| `--bam`       | `-b`       | String  | The path to the sorted BAM file (required).                                                   |
| `--output`    | `-o`       | String  | The path to the output file (required).                                                       |
| `--minreadlen`| `-m`       | Integer | Minimum read length to consider (optional). Default: 1000 for telomere capture, use 4000 for WGS. |
| `--maxgaplen` | `-g`       | Integer | Maximum allowed "gap" of non-canonical or 1bp mismatch sequence within telomere region (optional). Default: 250. |
| `--threads`   | `-t`       | Integer | Number of processing threads to use (optional). Default: Uses all available CPU threads.       |
| `--memlimit`   | `-l`       | Integer | Maximum RAM allocation per read batch (optional). Default: 8 Gb     |





