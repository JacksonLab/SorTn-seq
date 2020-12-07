# SorTn-seq: a functional genomics approach to discover regulators of bacterial gene expression

## Citation

Smith et al., 2020, Nature Microbiology: "The Rcs stress response inversely controls surface and CRISPR–Cas adaptive immunity to discriminate plasmids and phages"
Leah M. Smith, Simon A. Jackson, Lucia M. Malone, James E. Ussher, Paul P. Gardner and Peter C. Fineran*

Smith et al., unpublished: "SorTn-seq: a functional genomics approach to discover regulators of bacterial gene expression"
Leah M. Smith, Simon A. Jackson, Paul P. Gardner and Peter C. Fineran*

## SorTn-seq overview

SorTn-seq uses fluorescent reporters, saturation transposon mutagenesis and fluorescence activated cell sorting (FACS) to isolate bacterial mutants with altered gene expression. Sorted cell pools are deep sequenced to identify transposon insertion sites and the enrichment of mutants in high or low fluorescence bins is used to identify putative regulators of gene expression.

**This repository contains:**

1) The data analysis scripts from Smith et al., 2020, Nature Microbiology: SorTn-seq/Nature_Microbiology/

2) Data analysis scripts and an example dataset from the subsequent Protcol paper: Smith et al., unpublished.

### SorTn-seq data analysis overview



Figure 9. Summary of input and output files of the SorTn-seq analysis. FASTQ files are first processed to assess quality and remove adaptor contamination in the terminal window (shell). Processed files are fed into the TraDIS pipeline to identify the transposon tag and map reads to the reference genome (.fasta file). The TraDIS pipeline summarizes mapping and insertion statistics (.stats file), as well as producing sample-specific files, such reads per nucleotide position (.plot files) and Binary Alignment Map (BAM) files and indices (.bam and .bam.bai). In the terminal, BAM files are converted to Browser Extensible Data (BED) files (.bed) for subsequent analysis in R. To assign mapped reads to specific genomic features, an organism-specific feature table ([genome.prefix]_features_sortnseq.xlsx) is first generated in R (SorTnSeq_format_features.R), which parses RefSeq General Feature Format (GFF) files and adds intergenic regions as features. The feature table, BED files, and user-supplied sample information (sample_metadata.xlsx) are used to generate tables of read counts, insertion counts, and insertion index (number of insertions / feature length) for each sample (.xlsx files). To identify differentially enriched features, the unique insertion table (SorTnSeq_unique_insertions.xlsx) and insertion index table (SorTnSeq_table_insertion_index.xlsx) are processed using edgeR (SorTnSeq_analysis.R). A table summarizing feature enrichment (SorTnSeq_results_depleted_unique_insertions.xlsx) is generated along with plots that summarize the results (.pdf).


## Data analysis

**Process raw data:**
Requires the RefSeq nucleotide fasta file ([genome.prefix]_genomic.fna)

```bash
# Quality control of raw sequencing data
fastqc -t 32 *.fastq.gz

# Optional read trimming
trimmomatic SE -threads 20 -trimlog trim_summary [input].fastq.gz [output].fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:1

# Bio-TraDIS
find *.fastq.gz -printf '%f\n' > filelist.txt
bacteria_tradis --smalt --smalt_k 10 --smalt_s 1 --smalt_y 0.92 --smalt_r -1 -mm 2 -v -f filelist.txt -T TATAAGAGACAG -r [genome.prefix]_genomic.fna

# Convert .bam files to .bed format
for FILE in *.bam; do
bedtools bamtobed -i $FILE > $FILE.bed
done

```

**Generate a list of genome features and add intergenic regions.**
Requires the RefSeq .gff file () corresponding to the genome assembly used above.




## Dependencies:

- **FastQC**
*reference*
- **Trimmomatic**
*reference*
- **Bio-TraDIS**
*reference*
- **R**  
*Version 4.0.3 or higher https://www.r-project.org/*










