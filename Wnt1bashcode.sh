wget http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/Mus_musculus.GRCm39.105.gtf.gz
gtfToGenePred -genePredExt Mus_musculus.GRCm39.105.gtf /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > Mus_musculus.GRCm39.105.refflat

#indexing genome for HISAT2 alignment: indexing genome needs high computation.
hisat2-build -p 7 \
 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa \
 HISAT_index/Bulk_index

#steps from 13 to 25 should be run on all samples P14_K_I is just an example. (except FastQc and Multiqc all samples are included)   
 #map the data to the indexed genome 
fastqc -o FastQC *.fastq.gz
trim_galore --paired --nextera -o ../trimmed JP_46_S162_L004_R1_001.fastq.gz JP_46_S162_L004_R2_001.fastq.gz -j 7
fastqc -o FastQC *.fq.gz
 hisat2 -x ../../../../reference/Hisat_index/Bulk_index -1 JP_46_S162_L004_R1_001_val_1.fq.gz  -2 JP_46_S162_L004_R2_001_val_2.fq.gz -S ../../bam/P14_wnt1_I.sam -p 16 -t
 #convert SAM files to BAM files.
 samtools view -b -@ 8 P14_K_I.sam > P14_K_I.bam
 #sort BAM files.
 samtools sort -@ 8 P14_K_I.bam > P14_K_I.sorted.bam
 #index BAM files.
samtools index P14_K_I.sorted.bam
 #QC of the alignment; metrics with Picard Tools.
java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=P14_K_I.sorted.bam OUTPUT=P14_K_I.mkdup.bam METRICS_FILE=P14_K_I.mkdup_metrics.txt CREATE_INDEX=true
java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=P14_K_I.sorted.bam OUTPUT=P14_K_I.alignement_metrics.txt REFERENCE_SEQUENCE=../../../reference/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics INPUT=P14_Wnt1_I.sorted.bam OUTPUT=P14_Wnt1_I.insert_size.txt HISTOGRAM_FILE=P14_Wnt1_I.insert_size.pdf
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics INPUT=P14_Wnt1_I.sorted.bam REF_FLAT=../../../reference/Mus_musculus.GRCm39.105.refflat OUTPUT=P14_Wnt1_I.RNA_metrics.txt STRAND=NONE
 
 multiqc .
 
 #indexing transcriptome for Salmon. indexing transcriptome needs high computation
cat Mus_musculus.GRCm39.cdna.all.fa.gz Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz > gentrome.fa.gz
grep "^>" <(gunzip -c GRCm39.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
salmon index -t gentrome.fa.gz -d decoys.txt -k 15 -p 50 -i Salmon_index1 --gencode

 #Salmon Gene expression quantification :should be run on all samples.
  salmon quant -p 20 -i ../../../../reference/index_Salmon/Salmon_index1 --gcBias -l A -1 ../data/trimmed/JP_46_S162_L004_R1_001_val_1.fq.gz -2 ../data/trimmed/JP_46_S162_L004_R2_001_val_2.fq.gz -o ../Salmon_O/P14_Wnt1_I
 #Make transcript to gene table.
 zcat Mus_musculus.GRCm39.cdna.all.fa.gz | grep "^>" | cut -f 1,4 | sed -e 's/^>//' -e 's/gene://' -e 's/\.[0-9]*$//' | tr ' ' '\t' > tx2gene.tsv