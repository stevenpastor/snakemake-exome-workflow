# Simple Workflow on exome sequencing data.
## Paper: Clark M.J., et al. Performance comparison of exome DNA sequencing technologies. Nature biotechnology 2011; 29(10):908-914

### Obtain data. 
```
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar zxvf sratoolkit.current-centos_linux64.tar.gz
rm sratoolkit.current-centos_linux64.tar.gz

sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --split-3 SRR309293 # 34 total gigabytes.
sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --split-3 SRR309292 # this one is pretty large.
sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --split-3 SRR309291

# Reference from UCSC; be certain to read the text in the link for an understanding of the various file types:
http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/
```
* For simplicity, I did not include the unlocalized (chrUn) or unplaced (chr_random) sequences, only the chrNUM.fa.gz files.


### Use a submitter file to submit to compute node on Proteus.
* The file includes:
* Index Reference Genome, hg38, for faster searching with BWA alignment:
* command: bwa index <reference.fa>
* can use compressed (gzip) reference fasta file.
```
bwa index hg38.fa.gz
```

### Note how the reference is the "base" filename.
* no need to include the bwa indexed files in the alignment call.
* be certain the hg38.fa.gz is in same path as indexed files and bwa will automatically find indices.

### Map paired-end reads to indexed reference:
* PE command: bwa mem [options] <reference.fa> <fastq_mate_1.fq> <fastq_mate_2.fq> > <output.sam>
* I am using 64 threads here (-t 64) and in certain downstream commands (-@ 64 in Samtools).
```
bwa mem -t 64 hg38.fa.gz SRR309293_1.fastq SRR309293_2.fastq > SRR309293.sam
```

### Convert SAM to binary BAM, for speed of searching and space-saving:
```
samtools view -@ 64 -b SRR309293.sam > SRR309293.bam
```

### Sort the BAM for again, speed of searching:
```
samtools sort -@ 64 SRR309293.bam > SRR309293_sorted.bam
```

### Index sorted BAM for speed of searching:
```
samtools index SRR309293_sorted.bam
```

* no need to provide output for the index, as provided automatically (.bam.bai).

### To find reads mapped to particular region in the genome:
```
samtools view SRR309293_sorted.bam chr1:1000000-1000500 | head
```

* note: this prints first 10 lines only ("head").

### Summary of commands, with use of lustre and alignment in one line:
```
bwa index hg38.fa.gz
cp hg38.fa.gz* /lustre/scratch/<ourGrp>/<dir_of_choice>
cp SRR309293*fastq /lustre/scratch/<ourGrp>/<dir_of_choice>
cd /lustre/scratch/<ourGrp>/<dir_of_choice>

bwa mem -t 64 hg38.fa.gz SRR309293_1.fastq SRR309293_2.fastq | \
samtools view -@ 64 -b - | \
samtools sort -@ 64 - > SRR309293_sorted.bam

samtools index SRR309293_sorted.bam

cp SRR309293_sorted.bam <dir_of_choice>
```

* "|" is a pipe.
* Takes output of left command (left of "|") and uses as input in command to right of the "|".
* The "-" replaces the input file (for those familiar, this is also called standard in or stdin).
* Samtools recognizes "-" as stdin and stdout (where applicable).

### Call variants using samtools + bcftools:
```
samtools mpileup -g -f hg38.fa.gz SRR309293_sorted.bam | bcftools call -mv - > SRR309293.vcf
```

### (Optional) Convert to bed:
```
vcf2bed --snvs < SRR309293.vcf > SRR309293.bed
```
