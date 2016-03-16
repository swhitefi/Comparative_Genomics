# Bacterial Comparative Genomics Workshop

- [Day 1 Morning](https://github.com/alipirani88/Comparative_Genomics#day-1-morning)
	- [Getting your data onto Flux and setting up Environment variable](https://github.com/alipirani88/Comparative_Genomics#getting-your-data-onto-glux-and-setting-up-environment-variable)
	- [Quality Control using FastQC](https://github.com/alipirani88/Comparative_Genomics#quality-control-using-fastqc)
	- [Quality Trimming using Trimmomatic](https://github.com/alipirani88/Comparative_Genomics#quality-trimming-using-trimmomatic)

- [Day 1 Afternoon](https://github.com/alipirani88/Comparative_Genomics#day-1-afternoon)
	- [Read Mapping](https://github.com/alipirani88/Comparative_Genomics#read-mapping)
	- [Variant Calling](https://github.com/alipirani88/Comparative_Genomics#variant-calling-and-filteration)

# Day 1 Morning
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

## Getting your data onto Flux and setting up Environment variable

**Log in to Flux**

```
ssh user@flux-login.engin.umich.edu
```

> *user = your umich unique name

**Set up your .bashrc file so your environment is all set for genomic analysis!**

>i. Make a backup copy of ~/.bashrc file in case something goes wrong 
	
```
cp ~/.bashrc ./bashrc_backup
```
	
>ii. Add a line to your .bashrc file that points to required Perl library directories

```
export PERL5LIB=/scratch/micro612w16_fluxod/shared/bin/PAGIT/lib:/scratch/micro612w16_fluxod/shared/bin/vcftools_0.1.12b/perl:$PERL5LIB
```

>iii. Add entries in your .bashrc file to add genomics programs to your path variable

```
export PATH=$PATH: /scratch/micro612w16_fluxod/shared/bin/mauve_snapshot_2015-02-13/linux-x64/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/mauve_snapshot_2015-02-13/linux-x64/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/blast/bin/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/vcftools_0.1.12b/perl/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/tabix-0.2.6/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/bwa-0.7.12/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/bcftools-1.2/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/samtools-1.2/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/sratoolkit/bin/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/Spades/bin/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/FastQC/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/GenomeAnalysisTK-3.3-0/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/picard-tools-1.130/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/qualimap_v2.1/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/vcftools_0.1.12b/bin/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/sratoolkit/bin/
```

>iv. Source your .bashrc file

```
source .bashrc
```

## Quality Control using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ "FastQC homepage")
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

>i. Create directory for analysis

```
mkdir Rush_KPC_264_FastQC_results
mkdir Rush_KPC_264_FastQC_results/before_trimmomatic
```

>ii. Verify that FastQC is in your path

```
fastqc -h
```

>iii. Copy example fastq files to your home directory

```
cp /scratch/micro612w16_fluxod/shared/example_sample_2008_data_for_QC/Rush_KPC_264_*.gz ./
```

>iv. Get an interactive cluster node to start running programs

>v. Run FastQC to generate quality report on one of the example files

```
fastqc -o Rush_KPC_264_FastQC_results/before_trimmomatic/ Rush_KPC_264_1_combine.fastq.gz Rush_KPC_264_2_combine.fastq.gz –extract
```

>vi. Exit your cluster node so you don’t waste cluster resources and $$$!

>vii. Download FastQC report to your home computer to examine

`Screenshots explanation.`

>viii. Write shell script to run FastQC on all files in ???? directory(Not necessary)

>ix. Write PBS script to run your FastQC shell script on the cluster(Not necessary)

>x. Once your cluster job is finished, download reports to your home computer(Not necessary)

## Quality Trimming using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic "Trimmomatic Homepage")
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

>i. Get an interactive cluster node to start running programs


>ii. Create output directories to save trimmomatic results

```
mkdir Rush_KPC_264_trimmomatic_results
mkdir Rush_KPC_264_trimmomatic_results_with_headcrop/
mkdir Rush_KPC_264_FastQC_results/after_trimmomatic
mkdir Rush_KPC_264_FastQC_results/after_trimmomatic_headcrop/
```

>iii. Load latest version of java and try to run trimmomatic

```
module load lsa java/1.8.0

java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar –h
```

>-- explaining parameters and its default value. Adapter file. Changing only SLIDINGWINDOW parameter from default 4:15 to 4:20 for raw reads.

>iv. Run trimmomatic on raw reads

```
time java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_264_1_combine.fastq.gz Rush_KPC_264_2_combine.fastq.gz Rush_KPC_264_trimmomatic_results/forward_paired.fq.gz Rush_KPC_264_trimmomatic_results/forward_unpaired.fq.gz Rush_KPC_264_trimmomatic_results/reverse_paired.fq.gz Rush_KPC_264_trimmomatic_results/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:40 HEADCROP:0
```

>v. Run FastQC on trimmomatic results and check report on your local computer

```
fastqc -o Rush_KPC_264_FastQC_results/after_trimmomatic/ --extract -f fastq Rush_KPC_264_trimmomatic_results/forward_paired.fq.gz Rush_KPC_264_trimmomatic_results/reverse_paired.fq.gz
```
	
>-- explain fastqc results with screenshots. How head bases in per base sequence content graph are imbalanced? The cross signal sign for that graph? How you can fix it by using headcrop parameter in trimmomatic? Fastqc report.

>vi. Run trimmomatic with headcrop 9

```
time java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_264_1_combine.fastq.gz Rush_KPC_264_2_combine.fastq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/forward_paired.fq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/forward_unpaired.fq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/reverse_paired.fq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:40 HEADCROP:9
```
		
>-- explain per base sequence content changed to just warning from cross sign.

>vii. Run FastQC on updated trimmomatic results and check report on your local computer

```
fastqc -o Rush_KPC_264_FastQC_results/after_trimmomatic_headcrop/ --extract -f fastq Rush_KPC_264_trimmomatic_results_with_headcrop/forward_paired.fq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/reverse_paired.fq.gz
```

[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

# Day 1 Afternoon

## Read Mapping
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

**1. Create a directory to save results and run trimmomatic**

```
 mkdir Rush_KPC_264_varcall_result
```

```
 java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_264_1_combine.fastq.gz
Rush_KPC_264_2_combine.fastq.gz Rush_KPC_264_varcall_result/forward_paired.fq.gz Rush_KPC_264_varcall_result/forward_unpaired.fq.gz Rush_KPC_264_varcall_result/reverse_paired.fq.gz Rush_KPC_264_varcall_result/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:40
HEADCROP:0
```

**2. Map your reads against a finished reference genome using [BWA](http://bio-bwa.sourceforge.net/bwa.shtml "BWA manual")**

>i. Create BWA index from Reference fasta file:

```
bwa index /scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta
```

>index file usage 

>ii. Align reads to reference and output into SAM file

```
bwa mem -M -R "@RG\tID:96\tSM:Rush_KPC_264_1_combine.fastq.gz\tLB:1\tPL:Illumina" -t 8
/scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta Rush_KPC_264_varcall_result/forward_paired.fq.gz Rush_KPC_264_varcall_result/reverse_paired.fq.gz > Rush_KPC_264_varcall_result/Rush_KPC_264__aln.sam
```

> -R readgroup parameter; what does it say?

**3. SAM/BAM manipulation and variant calling using [Samtools](http://www.htslib.org/doc/samtools.html "Samtools Manual")

>i. Change directory to results folder:

```
cd Rush_KPC_264_varcall_result
```

>ii. Convert SAM to BAM using SAMTOOLS:

```
samtools view -Sb Rush_KPC_264__aln.sam > Rush_KPC_264__aln.bam
```

>iii. Sort BAM file using SAMTOOLS:

```
samtools sort Rush_KPC_264__aln.bam Rush_KPC_264__aln_sort
```

**4. Mark duplicates(PCR optical duplicates) and remove them using [PICARD](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates "Picard MarkDuplicates")**

>i. Create a dictionary for reference fasta file required by PICARD(If KPNIH1.dict doesn’t exist).
 
```
java -jar /scratch/micro612w16_fluxod/shared/bin/picard-tools-1.130/picard.jar CreateSequenceDictionary REFERENCE=/scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta OUTPUT=/scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.dict
```

>ii. Run PICARD for removing duplicates.

```
java -jar /scratch/micro612w16_fluxod/shared/bin/picard-tools-1.130/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=Rush_KPC_264__aln_sort.bam OUTPUT= Rush_KPC_264__aln_marked.bam METRICS_FILE=Rush_KPC_264__markduplicates_metrics
CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
```

>PCR optical duplicates; why remove that?

>iii. Sort these marked BAM file again (for downstream compatibility)

```
samtools sort Rush_KPC_264__aln_marked.bam Rush_KPC_264__aln_sort
```

>iv. Index these marked bam file using SAMTOOLS

```
samtools index Rush_KPC_264__aln_sort.bam
```

## Variant Calling and Filteration
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

**1. Call variants using [samtools](http://www.htslib.org/doc/samtools.html "samtools manual") mpileup and [bcftools](https://samtools.github.io/bcftools/bcftools.html "bcftools")**

```
samtools mpileup -ug -f /scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta Rush_KPC_264__aln_sort.bam | bcftools call -O v -v -c -o Rush_KPC_264__aln_mpileup_raw.vcf
```

>**-g generate genotype likelihood in bcf format   
>**mpileup format   
>**-c samtools consensus caller

**2. Variant filtering and processed file generation using GATK and vcftools**

>i. Variant filtering using [GATK](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php "GATK Variant Filteration"):
 
```
java -jar /scratch/micro612w16_fluxod/shared/bin/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T VariantFiltration -R
/home2/apirani/bin/reference/KPNIH1/KPNIH1.fasta -o Rush_KPC_264__filter_gatk.vcf --variant Rush_KPC_264__aln_mpileup_raw.vcf --filterExpression "FQ < 0.025 && MQ > 50 && QUAL > 100 && DP > 15" --filterName pass_filter
```

>**FQ, MQ, DP, QUAL

>ii. Remove indels and keep only variants that passed filter parameter from VCF file using [vcftools](http://vcftools.sourceforge.net/man_latest.html vcftools manual):
 
 
```
vcftools --vcf Rush_KPC_264__filter_gatk.vcf --keep-filtered pass_filter --remove-indels --recode --recode-INFO-all --out
Rush_KPC_264__filter_onlysnp
```

>>**why remove indels  

>iii. Generate Consensus fasta file from filtered variants using vcftools:

```
bgzip Rush_KPC_264__filter_onlysnp.recode.vcf
tabix Rush_KPC_264__filter_onlysnp.recode.vcf.gz
```










