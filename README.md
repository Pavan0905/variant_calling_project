# Variant Calling and Annotation Pipeline

## Overview
This project processes raw sequencing data to identify and annotate high-impact genetic variants, particularly in the **ARID1A** gene. The pipeline leverages industry-standard bioinformatics tools for alignment, variant calling, filtering, annotation, and visualization.

## Directory Structure
```
variant_calling_project/
│-- data/                      # Contains raw and processed sequencing data
│-- results/                   # Stores output files, figures, and final reports
│   ├── fastqc_reports/        # Quality control reports
│   ├── variant_calling/       # VCF files and variant annotations
│   ├── figures/               # Visualizations
│-- scripts/                   # Analysis scripts
│-- docs/                      # Documentation and project report
│-- README.md                  # Project overview and instructions
```

## Tools & Technologies Used
- **FastQC** - Quality control of raw sequencing reads
- **Trimmomatic** - Adapter and quality trimming
- **BWA-MEM** - Read alignment to the reference genome
- **Samtools** - BAM file processing and variant filtering
- **GATK** - Variant calling and filtering
- **SnpEff** - Variant annotation
- **Python (Matplotlib, Pandas)** - Data analysis and visualization
- **GitHub** - Version control and project hosting

## Installation
Ensure the required tools are installed using:
```bash
sudo apt-get install fastqc bwa samtools
conda install -c bioconda gatk4 snpeff trimmomatic
```

## Pipeline Workflow

### **Step 1: Quality Control**
Run FastQC to assess sequencing read quality:
```bash
fastqc data/raw_reads/*.fastq -o results/fastqc_reports/
```
Trim adapters and low-quality reads using Trimmomatic:
```bash
trimmomatic PE -threads 4 \
  data/raw_reads/SRR925811_1.fastq data/raw_reads/SRR925811_2.fastq \
  data/processed_reads/SRR925811_1_paired.fastq data/processed_reads/SRR925811_1_unpaired.fastq \
  data/processed_reads/SRR925811_2_paired.fastq data/processed_reads/SRR925811_2_unpaired.fastq \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
```

### **Step 2: Read Alignment**
Align reads to the reference genome using BWA-MEM:
```bash
bwa mem -t 4 data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  data/processed_reads/SRR925811_1_paired.fastq data/processed_reads/SRR925811_2_paired.fastq > data/aligned_reads/SRR925811_aligned.sam
```
Convert to BAM, sort, and index:
```bash
samtools view -bS data/aligned_reads/SRR925811_aligned.sam > data/aligned_reads/SRR925811_aligned.bam
samtools sort data/aligned_reads/SRR925811_aligned.bam -o data/aligned_reads/SRR925811_sorted.bam
samtools index data/aligned_reads/SRR925811_sorted.bam
```

### **Step 3: Variant Calling**
Use GATK HaplotypeCaller to generate a VCF file:
```bash
gatk HaplotypeCaller -R data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -I data/aligned_reads/SRR925811_sorted.bam -O results/variant_calling/SRR925811_variants.g.vcf
```

### **Step 4: Variant Filtering and Annotation**
Filter **high-confidence** variants using GATK SelectVariants:
```bash
gatk SelectVariants -R data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -V results/variant_calling/SRR925811_variants.g.vcf --select-type-to-include SNP \
  -O results/variant_calling/SRR925811_filtered_variants.vcf
```
Annotate variants using **SnpEff**:
```bash
java -jar snpEff.jar GRCh38.86 results/variant_calling/SRR925811_filtered_variants.vcf \
  > results/variant_calling/SRR925811_annotated_variants.vcf
```

### **Step 5: Visualization & Analysis**
Analyze and visualize **high-impact genes** using Python:
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/variant_calling/SRR925811_annotated_variants.vcf", sep='\t', comment='#', header=None)
gene_counts = df[7].str.extract(r'ANN=.*?\|HIGH\|.*?\|.*?\|(.*?)\|')[0].value_counts().head(10)

plt.figure(figsize=(10,5))
gene_counts.plot(kind='bar', color='steelblue')
plt.xlabel("Gene")
plt.ylabel("Number of High-Impact Variants")
plt.title("Top Genes with High-Impact Variants")
plt.xticks(rotation=45)
plt.savefig("results/figures/high_impact_genes_chart.png")
plt.show()
```

## Results
- **Total Variants Identified:** **124,117**
  - **SNVs:** 114,962
  - **Insertions:** 4,178
  - **Deletions:** 4,971
- **Key Genes with High-Impact Variants:** **ARID1A, CELA1, OR8U1, MIA2, FCGR2C**

## Repository & Reproducibility
This project is hosted on **GitHub**. Since VCF files are large, only placeholder files are included.
To retrieve full datasets:
```bash
wget [Data_Source_URL]
```
To run the analysis:
```bash
bash scripts/run_pipeline.sh
```

## Future Improvements
- Integrate Nextflow/Snakemake for automation
- Improve annotation using additional databases (ClinVar, COSMIC)
- Deploy results using **R Shiny** or **Python Streamlit**

## Contributors
- **Pavan Kumar Singuru**  
  Northeastern University, US

## Contact
For questions or collaboration, reach out via **GitHub Issues** or email at singuru.p@domain.com

