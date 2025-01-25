# Variant Calling Project

## Objective
Analysis of high-impact variants in the ARID1A gene and other top genes, using variant annotation and visualization.


## Data Used

### 1. Annotated Variants File
- **File Name**: `SRR925811_annotated_variants.vcf`
- **Source**: Publicly available dataset from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra), accession number `SRR925811`.
- **Description**: Annotated variant calls, including SNVs, insertions, and deletions, across the human genome.

### 2. Raw Reads
- **Files**:
  - `data/raw_reads/SRR925811_1.fastq`
  - `data/raw_reads/SRR925811_2.fastq`
- **Source**: Same NCBI SRA dataset as above.
- **Description**: Paired-end raw reads.

### 3. Reference Genome
- **File Name**: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- **Source**: [Ensembl Genome Browser](https://www.ensembl.org/).
- **Description**: Human genome assembly GRCh38 used for alignment and variant calling.



## Steps Performed
1. **Data Organization**:
   - Created a structured project directory with `data/`, `results/`, `scripts/`, and `docs/` folders.
   - Stored the annotated VCF file (`SRR925811_annotated_variants.vcf`) in the `data/` folder.
   - Stored analysis scripts (e.g., `variant_analysis.py`) in the `scripts/` folder.
   - Saved results, including figures, in the `results/` folder.


## Workflow
1. **Data Preprocessing**
   - Performed quality checks using FastQC.
   - Trimmed adapters and low-quality bases using Trimmomatic.
2. **Read Alignment**
   - Aligned reads to the human reference genome (GRCh38) using BWA.
3. **Variant Calling and Annotation**
   - Detected variants using GATK.
   - Filtered variants and identified high-impact variants.
4. **Visualization**
   - Created visualizations of high-impact genes using Matplotlib.






2. **Variant Analysis**:
   - Extracted and counted the number of high-impact variants specific to the ARID1A gene.
   - Categorized ARID1A variants by type:
     - **SNVs**: 2
     - **Insertions (INS)**: 0
     - **Deletions (DEL)**: 0
   - Identified and summarized high-impact genes across the dataset.

3. **Visualization**:
   - Generated a bar chart (`high_impact_genes_chart.png`) showing the number of high-impact variants for the top genes, including ARID1A.
   - Visualized ARID1A-specific variant types for better understanding.

4. **Summary Statistics**:
   - Total variants analyzed: **124,117**.
   - Breakdown by type:
     - **SNVs**: 114,962
     - **Insertions (INS)**: 4,178
     - **Deletions (DEL)**: 4,971
     - **Others**: 6
   - High-impact genes identified: **ARID1A**, **CELA1**, **OR8U1**, and others.

## Results
- **ARID1A Variant Statistics**:
  - SNVs: 2
  - Insertions (INS): 0
  - Deletions (DEL): 0
- Generated the `high_impact_genes_chart.png`, visualizing the top genes with high-impact variants.
- Highlighted genes with the highest number of high-impact variants:
  - **CELA1**, **OR8U1**, **MIA2**, **FCGR2C**, and **ARID1A**.

## Next Steps
1. Perform functional annotation of high-impact variants to assess their potential effects on gene function and disease association.
2. Extend the analysis to additional genes or regions of interest.
3. Automate the pipeline using workflow tools like Nextflow or Snakemake for scalability.
4. Integrate the results with external variant databases (e.g., ClinVar, COSMIC) to enhance biological interpretation.
5. Publish findings in a report or dashboard using visualization tools like R Shiny or Python Streamlit.


## Resources
- [Bash Commands Log](bash_commands_log.md)
- [Tools Used](tools_used.md)
