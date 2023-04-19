#!/bin/bash
trimmomatic PE -threads 8 -phred33 /home/ubuntu/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/untrimmed_fastq/NGS0001.R1.fastq.gz /home/ubuntu/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/untrimmed_fastq/NGS0001.R2.fastq.gz -baseout /home/ubuntu/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/trimmed_fastq/NGS0001_trimmed_R.fastq ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:15

cd /home/ubuntu/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/trimmed_fastq/

fastqc -t 8 *.fastq

mkdir ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/fastqc_trimmed_reads
 
mv *fastqc* ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/fastqc_trimmed_reads
 
cd ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/fastqc_trimmed_reads
 
for zip in *.zip; do unzip $zip; done
 
head NGS0001_trimmed_R_1P_fastqc//summary.txt

head NGS0001_trimmed_R_1U_fastqc/summary.txt

 
head NGS0001_trimmed_R_2P_fastqc/summary.txt

 
head NGS0001_trimmed_R_2U_fastqc/summary.txt

 
cat */summary.txt > ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSlogs/fastqc_trimmed_summaries.txt
 

mkdir -p ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference
 
mv ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/hg19.fa.gz ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/

bwa index ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/hg19.fa.gz

ls ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/


cd ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/

chmod +x hg19.fa.gz

ls

mkdir ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/aligned_data

bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera_NGS0001_blood\tDT:2023-04-12\tPU:11V6WR1' -I 250,50  ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/hg19.fa.gz ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/fastqc_trimmed_reads/NGS0001_trimmed_R_1P.fastq ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/fastqc_trimmed_reads/NGS0001_trimmed_R_2P.fastq > ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/aligned_data/NGS0001.bam

samtools sort NGS0001.bam > NGS0001_sorted.bam

samtools index NGS0001_sorted.bam


mv /home/ubuntu/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/aligned_data/NGS0001_sorted.bam.bai ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/

cd Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/aligned_data/

conda innit new-env
conda activate new-env

conda install -c bioconda picard

picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

conda deactivate

samtools index NGS0001_sorted_marked.bam


samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam

samtools index NGS0001_sorted_filtered.bam

samtools flagstat NGS0001_sorted_marked.bam >> NGS0001_sorted_marked_flagstat.txt

samtools flagstat NGS0001_sorted_filtered.bam >> NGS0001_sorted_filtered_flagstat.txt

samtools view NGS0001_sorted.bam | less -S

samtools view NGS0001_sorted_filtered.bam | less -S

samtools view NGS0001_sorted_marked.bam | less -S

samtools idxstats NGS0001_sorted_filtered.bam >> NGS0001_sorted_filtered_idxstats.txt


conda activate new-env

picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5

conda deactivate

pwd

cd /home/ubuntu/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/

gunzip ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/hg19.fa.gz > ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/hg19.fa

samtools faidx ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/hg19.fa

freebayes --bam ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/reference/hg19.fa --vcf ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001.vcf

sudo apt install tabix


bgzip ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001.vcf

tabix -p vcf ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001.vcf.gz


conda activate new-env


vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \                          
> ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001.vcf.gz > ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.vcf


bedtools intersect -header -wa -a ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.vcf -b ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSdata/annotation.bed ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.vcf

tabix -p vcf ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.vcf.gz


conda deactivate

cd annovar

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
 
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
 
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
 
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
 
./convert2annovar.pl -format vcf4 ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.vcf.gz > ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.avinput

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/

 
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

 
./convert2annovar.pl -format vcf4 ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.vcf.gz > ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.avinput

chmod +x table_annovar.pl
chmod +x coding_change.pl


./table_annovar.pl ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.avinput ~/annovar/humandb/ -buildver hg19  -out ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered

cd snpEff/

java -Xmx8g -jar snpEff.jar -c ~/snpEff/snpEff.config -v GRCh37.75 ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/NGS0001_filtered.vcf.gz > ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/test.NGS0001.filtered.ann.vcf

java -jar SnpSift.jar filter "ANN[0].EFFECT has 'exon_variant'" ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/test.NGS0001.filtered.ann.vcf > ~/Adv_Bioinformatics_Assignment/DNAseqAnalysis/NGSresults/test.NGS0001.ann.filter_exon_first.vcf

