# Bash Commands Log

# Navigate to project directory
cd ~/variant_calling_project

# Create data folder and move VCF file
mkdir -p data
mv ~/Downloads/SRR925811_annotated_variants.vcf data/

# Count variants for ARID1A
grep "ARID1A" data/SRR925811_annotated_variants.vcf | \
awk -F'\t' '{if (length($4) == 1 && length($5) == 1) print "SNV"; else if (length($4) < length($5)) print "INS"; else if (length($4) > length($5)) print "DEL"; else print "OTHER"}' | \
sort | uniq -c

# Generate chart of high-impact genes
python scripts/variant_analysis.py





### Data Preprocessing
```bash
# FastQC quality check
fastqc -o results/fastqc_reports/ data/raw_reads/SRR925811_1.fastq data/raw_reads/SRR925811_2.fastq

# Adapter trimming with Trimmomatic
trimmomatic PE -phred33 \
  data/raw_reads/SRR925811_1.fastq data/raw_reads/SRR925811_2.fastq \
  data/processed_reads/SRR925811_1_paired.fastq data/processed_reads/SRR925811_1_unpaired.fastq \
  data/processed_reads/SRR925811_2_paired.fastq data/processed_reads/SRR925811_2_unpaired.fastq \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50



Read Alignment
# Align reads to the reference genome
bwa mem -t 4 data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  data/processed_reads/SRR925811_1_paired.fastq data/processed_reads/SRR925811_2_paired.fastq > data/aligned_reads/SRR925811_aligned.sam

# Convert SAM to BAM, sort, and index
samtools view -Sb data/aligned_reads/SRR925811_aligned.sam > data/aligned_reads/SRR925811_aligned.bam
samtools sort data/aligned_reads/SRR925811_aligned.bam > data/aligned_reads/SRR925811_sorted.bam
samtools index data/aligned_reads/SRR925811_sorted.bam


Variant Calling and Filtering
# Call variants with GATK
gatk HaplotypeCaller -R data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -I data/aligned_reads/SRR925811_sorted.bam \
  -O results/variant_calling/SRR925811_variants.g.vcf

# Filter high-quality variants
gatk VariantFiltration -R data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -V results/variant_calling/SRR925811_variants.g.vcf \
  -O results/variant_calling/SRR925811_filtered_variants.vcf

Visualization
# Generate a bar chart of high-impact genes
python scripts/variant_analysis.py

