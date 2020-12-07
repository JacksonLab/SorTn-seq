# SorTn-seq: a functional genomics approach to discover regulators of bacterial gene expression

## Citation

Data analysis scripts from Smith et al., 2020: "The Rcs stress response inversely controls surface and CRISPR–Cas adaptive immunity to discriminate plasmids and phages"
Leah M. Smith, Simon A. Jackson, Lucia M. Malone, James E. Ussher, Paul P. Gardner and Peter C. Fineran*

Data analysis scripts and example dataset from Smith et al., unpublished: "SorTn-seq: a functional genomics approach to discover regulators of bacterial gene expression"
Leah M. Smith, Simon A. Jackson, Paul P. Gardner and Peter C. Fineran*


## SorTn-seq overview










### Data analysis overview



Figure 9. Summary of input and output files of the SorTn-seq analysis. FASTQ files are first processed to assess quality and remove adaptor contamination in the terminal window (shell). Processed files are fed into the TraDIS pipeline to identify the transposon tag and map reads to the reference genome (.fasta file). The TraDIS pipeline summarizes mapping and insertion statistics (.stats file), as well as producing sample-specific files, such reads per nucleotide position (.plot files) and Binary Alignment Map (BAM) files and indices (.bam and .bam.bai). In the terminal, BAM files are converted to Browser Extensible Data (BED) files (.bed) for subsequent analysis in R. To assign mapped reads to specific genomic features, an organism-specific feature table ([genome.prefix]_features_sortnseq.xlsx) is first generated in R (SorTnSeq_format_features.R), which parses RefSeq General Feature Format (GFF) files and adds intergenic regions as features. The feature table, BED files, and user-supplied sample information (sample_metadata.xlsx) are used to generate tables of read counts, insertion counts, and insertion index (number of insertions / feature length) for each sample (.xlsx files). To identify differentially enriched features, the unique insertion table (SorTnSeq_unique_insertions.xlsx) and insertion index table (SorTnSeq_table_insertion_index.xlsx) are processed using edgeR (SorTnSeq_analysis.R). A table summarizing feature enrichment (SorTnSeq_results_depleted_unique_insertions.xlsx) is generated along with plots that summarize the results (.pdf).


## Data analysis

```bash
# Quality control of raw sequencing data
fastqc -t 32 *.fastq.gz
```



