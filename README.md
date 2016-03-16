# Comparative_Genomics

- [Day 1 Morning](https://github.com/alipirani88/Comparative_Genomics#Day-1-Morning)
	- [Getting your data onto Flux and setting up Environment variable](https://github.com/alipirani88/Comparative_Genomics#Getting-your-data-onto-Flux-and-setting-up-Environment-variable)
	- [Quality Control using FastQC](https://github.com/alipirani88/Comparative_Genomics#Quality-Control-using-FastQC)
	- [Quality Trimming using Trimmomatic](https://github.com/alipirani88/Comparative_Genomics#Quality-Trimming-using-Trimmomatic)


























# Day 1 Morning

##Getting your data onto Flux and setting up Environment variable

**Log in to Flux**

`ssh user@flux-login.engin.umich.edu`

> *user = your umich unique name

###Set up your .bashrc file so your environment is all set for genomic analysis!

>i. Make a backup copy of ~/.bashrc file in case something goes wrong 
	
`cp ~/.bashrc ./bashrc_backup`
	
>ii. Add a line to your .bashrc file that points to required Perl library directories

`export PERL5LIB=/scratch/micro612w16_fluxod/shared/bin/PAGIT/lib:/scratch/micro612w16_fluxod/shared/bin/vcftools_0.1.12b/perl:$PERL5LIB`

>iii. Add entries in your .bashrc file to add genomics programs to your path variable

`export PATH=$PATH: /scratch/micro612w16_fluxod/shared/bin/mauve_snapshot_2015-02-13/linux-x64/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/mauve_snapshot_2015-02-13/linux-x64/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/blast/bin/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/vcftools_0.1.12b/perl/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/tabix-0.2.6/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/bwa-0.7.12/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/bcftools-1.2/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/samtools-1.2/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/sratoolkit/bin/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/Spades/bin/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/FastQC/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/GenomeAnalysisTK-3.3-0/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/picard-tools-1.130/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/qualimap_v2.1/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/vcftools_0.1.12b/bin/`
`export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/sratoolkit/bin/`

>iv. Execute commands in your .bashrc file

`source .bashrc`

## Quality Control using FastQC

>i. Create directory for analysis

`mkdir Rush_KPC_264_FastQC_results`
`mkdir Rush_KPC_264_FastQC_results/before_trimmomatic`

>ii. Verify that FastQC is in your path

`fastqc -h`

>iii. Copy example fastq files to your home directory

`cp /scratch/micro612w16_fluxod/shared/example_sample_2008_data_for_QC/Rush_KPC_264_*.gz ./`

>iv. Get an interactive cluster node to start running programs

>v. Run FastQC to generate quality report on one of the example files

`fastqc -o Rush_KPC_264_FastQC_results/before_trimmomatic/ Rush_KPC_264_1_combine.fastq.gz Rush_KPC_264_2_combine.fastq.gz –extract`

>vi. Exit your cluster node so you don’t waste cluster resources and $$$!

>vii. Download FastQC report to your home computer to examine

`Screenshots explanation.`

>viii. Write shell script to run FastQC on all files in ???? directory(Not necessary)

>ix. Write PBS script to run your FastQC shell script on the cluster(Not necessary)

>x. Once your cluster job is finished, download reports to your home computer(Not necessary)

##Quality Trimming using Trimmomatic

>i. Get an interactive cluster node to start running programs


>ii. Create output directories to save trimmomatic results

`mkdir Rush_KPC_264_trimmomatic_results`
`mkdir Rush_KPC_264_trimmomatic_results_with_headcrop/`
`mkdir Rush_KPC_264_FastQC_results/after_trimmomatic`
`mkdir Rush_KPC_264_FastQC_results/after_trimmomatic_headcrop/`

>iii. Load latest version of java and try to run trimmomatic

`module load lsa java/1.8.0`

`java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar –h`

>-- explaining parameters and its default value. Adapter file. Changing only SLIDINGWINDOW parameter from default 4:15 to 4:20 for raw reads.

>iv. Run trimmomatic on raw reads

`time java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_264_1_combine.fastq.gz Rush_KPC_264_2_combine.fastq.gz Rush_KPC_264_trimmomatic_results/forward_paired.fq.gz Rush_KPC_264_trimmomatic_results/forward_unpaired.fq.gz Rush_KPC_264_trimmomatic_results/reverse_paired.fq.gz Rush_KPC_264_trimmomatic_results/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:40 HEADCROP:0`

>v. Run FastQC on trimmomatic results and check report on your local computer

`fastqc -o Rush_KPC_264_FastQC_results/after_trimmomatic/ --extract -f fastq Rush_KPC_264_trimmomatic_results/forward_paired.fq.gz Rush_KPC_264_trimmomatic_results/reverse_paired.fq.gz`
	
>-- explain fastqc results with screenshots. How head bases in per base sequence content graph are imbalanced? The cross signal sign for that graph? How you can fix it by using headcrop parameter in trimmomatic? Fastqc report.

>vi. Run trimmomatic with headcrop 9

`time java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_264_1_combine.fastq.gz Rush_KPC_264_2_combine.fastq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/forward_paired.fq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/forward_unpaired.fq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/reverse_paired.fq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:40 HEADCROP:9`
		
>-- explain per base sequence content changed to just warning from cross sign.

>vii. Run FastQC on updated trimmomatic results and check report on your local computer

`fastqc -o Rush_KPC_264_FastQC_results/after_trimmomatic_headcrop/ --extract -f fastq Rush_KPC_264_trimmomatic_results_with_headcrop/forward_paired.fq.gz Rush_KPC_264_trimmomatic_results_with_headcrop/reverse_paired.fq.gz`


