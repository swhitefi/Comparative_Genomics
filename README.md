# Bacterial Comparative Genomics Workshop

- [Day 1 Morning](https://github.com/alipirani88/Comparative_Genomics#day-1-morning)
***
	- [Getting your data onto Flux and setting up Environment variable](https://github.com/alipirani88/Comparative_Genomics#getting-your-data-onto-glux-and-setting-up-environment-variable)
	- [Quality Control using FastQC](https://github.com/alipirani88/Comparative_Genomics#quality-control-using-fastqc)
	- [Quality Trimming using Trimmomatic](https://github.com/alipirani88/Comparative_Genomics#quality-trimming-using-trimmomatic)

- [Day 1 Afternoon](https://github.com/alipirani88/Comparative_Genomics#day-1-afternoon)
***
	- [Read Mapping](https://github.com/alipirani88/Comparative_Genomics#read-mapping)
	- [Variant Calling](https://github.com/alipirani88/Comparative_Genomics#variant-calling-and-filteration)
	- [Visualize BAM/VCF files in IGV/ACT](https://github.com/alipirani88/Comparative_Genomics#visualize-bam-and-vcf-files-in-igv-or-act)

- [Day 2 Morning](https://github.com/alipirani88/Comparative_Genomics#day-2-morning)
***
	- [Genome Assembly](https://github.com/alipirani88/Comparative_Genomics#genome-assembly)
	- [Assembly evaluation](https://github.com/alipirani88/Comparative_Genomics#assembly-evaluation)
	- [Compare assembly to reference genome and Post-assembly genome improvement](https://github.com/alipirani88/Comparative_Genomics#compare-assembly-to-reference-genome-and-post-assembly-genome-improvement)
	- [Map reads to the final ordered assembly](https://github.com/alipirani88/Comparative_Genomics#map-reads-to-the-final-ordered-assembly)
	- [Genome Annotation](https://github.com/alipirani88/Comparative_Genomics#genome-annotation)
	- [Visualize multiple assemblies](https://github.com/alipirani88/Comparative_Genomics#visualize-multiple-assemblies)

- [Day 2 Afternoon](https://github.com/alipirani88/Comparative_Genomics#day-2-afternoon)
***
	- [Determine which genomes contain beta-lactamase genes](https://github.com/alipirani88/Comparative_Genomics#determine-which-genomes-contain-beta-lactamase-genes)
	- [Identification of antibiotic resistance genes with LS-BSR and the ARDB database](https://github.com/alipirani88/Comparative_Genomics#identification-of-antibiotic-resistance-genes-with-ls-bsr-and-the-ardb-database)
	- [Perform pan-genome analysis with LS-BSR](https://github.com/alipirani88/Comparative_Genomics#perform-pan-genome-analysis-with-ls-bsr)
	- [Perform genome comparisons with ACT](https://github.com/alipirani88/Comparative_Genomics#perform-genome-comparisons-with-act)

- [Day 3 Morning](https://github.com/alipirani88/Comparative_Genomics#day-3-morning)
***
	- [Perform whole genome alignment with Mauve](https://github.com/alipirani88/Comparative_Genomics#perform-whole-genome-alignment-with-Mauve)
	- [Perform DNA sequence comparisons and phylogenetic analysis in ape(an R package)](https://github.com/alipirani88/Comparative_Genomics#perform-dna-sequence-comparisons-and-phylogenetic-analysis-in-ape)
	- [Perform SNP density analysis to discern evidence of recombination](https://github.com/alipirani88/Comparative_Genomics#perform-snp-density-analysis-to-discern-evidence-of-recombination)
	- [Perform recombination filtering with gubbins](https://github.com/alipirani88/Comparative_Genomics#perform-recombination-filtering-with-gubbins)
	- [Create annotated publication quality trees with iTOL](https://github.com/alipirani88/Comparative_Genomics#create-annotated-publication-quality-trees-with-itol)

- [Day 3 Afternoon](https://github.com/alipirani88/Comparative_Genomics#day-3-afternoon)
***
	- [Perform QC on fastq files](https://github.com/alipirani88/Comparative_Genomics#perform-qc-on-fastq-files)
	- [Examine results of SPANDx pipeline](https://github.com/alipirani88/Comparative_Genomics#examine-results-of-spandx-pipeline)
	- [Recombination detection and tree generation](https://github.com/alipirani88/Comparative_Genomics#recombination-detection-and-tree-generation)
	- [Phylogenetic tree annotation and visualization](https://github.com/alipirani88/Comparative_Genomics#phylogenetic-tree-annotation-and-visualization)
	- [Assessment of genomic deletions](https://github.com/alipirani88/Comparative_Genomics#assessment-of-genomic-deletions)

- [Helpful resources for microbial genomics](https://github.com/alipirani88/Comparative_Genomics#helpful-resources-for-microbial-genomics)
***






# Day 1 Morning
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

## Getting your data onto Flux and setting up Environment variable

**Log in to Flux**

```
ssh user@flux-login.engin.umich.edu
```

> *user = your umich unique name

**Set up your .bashrc file so your environment is all set for genomic analysis!**

Environment variables are a way of passing information from the shell to programs when you run them. Programs look "in the environment" for particular variables and if they are found will use the values stored. Some are set by the system, others by you, yet others by the shell, or any program that loads another program. All the softwares/tools that we need in this workshop are installed in a directory "/scratch/micro612w16_fluxod/shared/bin/". We will set the environment variable PATH in .bashrc file by exporting the required paths.

>i. Make a backup copy of ~/.bashrc file in case something goes wrong. 
	
```
cp ~/.bashrc ./bashrc_backup
```
	
>ii. Add a line to your .bashrc file that points to required Perl library directories.

```
export PERL5LIB=/scratch/micro612w16_fluxod/shared/bin/PAGIT/lib:/scratch/micro612w16_fluxod/shared/bin/vcftools_0.1.12b/perl:$PERL5LIB
```

>iii. Add entries in your .bashrc file to add required genomics programs to your path variable.

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
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/snpEff/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/PAGIT/ABACAS/
export PATH=$PATH:/scratch/micro612w16_fluxod/shared/bin/blast-2.2.26/bin/
```

>iv. Source your .bashrc file

```
source .bashrc
```

## Quality Control using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ "FastQC homepage")
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

>i. Execute the following commands to copy files for this morning’s exercises to your scratch directory: 

```
cd /scratch/micro612w16_fluxod/username
cp -r /scratch/micro612w16_fluxod/shared/data/day1_morn/ ./
cd /scratch/micro612w16_fluxod/username/day1_morn/
ls
```

As soon as you receive your sample data from sequencing centre, the first thing you do is check its quality using quality control tool such as FastQC. But before carrying out extensive QC, you can run a bash one-liner to get some basic statistics about the raw reads.

Run the following command to print total number of reads in each file, total number of unique reads, percentage of unique reads, most abundant sequence(useful to find adapter sequences or contamination), its frequency, and frequency of that sequence as a proportion of the total reads.

```
for i in *.gz; do zcat $i | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'; done
```

You can find more of such bash one-liners at Stephen Turner's github [page.](https://github.com/stephenturner/oneliners)

Now we will run FastQC on these raw data to assess its quality. FastQC is a quality control tool that reads in sequence data in a variety of formats(fastq, bam, sam) and can either provide an interactive application to review the results or create an HTML based report which can be integrated into any pipeline. It is generally the first step that you take upon receiving the sequence data from sequencing facility to get a quick sense of its quality or whether it exhibits any unusual properties(contamination or interesting biological features)

>ii. In your day1_morn directory, create a new directory for saving FastQC results.

```
mkdir Rush_KPC_266_FastQC_results
mkdir Rush_KPC_266_FastQC_results/before_trimmomatic
```

>iii. Verify that FastQC is in your path by invoking it from command line.

```
fastqc -h
```

FastQC can be run in two modes: "command line" or as a GUI (graphical user interface). We will be using command line version of it.

>iv. Get an interactive cluster node to start running programs

>v. Run FastQC to generate quality report of sequence reads.

```
fastqc -o Rush_KPC_266_FastQC_results/before_trimmomatic/ Rush_KPC_266_1_combine.fastq.gz Rush_KPC_266_2_combine.fastq.gz --extract
```

This will generate two results directory, Rush_KPC_266_1_combine_fastqc and Rush_KPC_266_2_combine_fastqc in output folder provided with -o flag. 
The summary.txt file in these directories indicates if the data passed different quality control tests. 
You can visualize and assess the quality of data by opening html report in a local browser.


>vi. Exit your cluster node so you don’t waste cluster resources and $$$!

>vii. Download FastQC report to your home computer to examine

```
scp username@flux-xfer.engin.umich.edu:/scratch/micro612w16_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/before_trimmomatic/*.html /path-to-local-directory/
```

The analysis in FastQC is performed by a series of analysis modules. The left hand side of the main interactive display or the top of the HTML report show a summary of the modules which were run, and a quick evaluation of whether the results of the module seem entirely normal (green tick), slightly abnormal (orange triangle) or very unusual (red cross). 

`Screenshots explanation.`
`Explaining Summary results, Basic statistics, per base sequence quality, overrepresented sequences(adapters) from before trimmomatic report.`

Notice the quality drop(per base sequence quality graph) at the end of Rush_KPC_266_2_combine_fastqc.html report. This is commonly observed in illumina samples. The reason for this drop is that as the number of sequencing cycles performed increases, the average quality of the base calls, as reported by the Phred Scores produced by the sequencer falls. 

Now, Check the overrepresented sequences graph and the kind of adapters that were used for sequencing these samples.(Truseq or Nextera)

Check out [this](https://sequencing.qcfail.com/articles/loss-of-base-call-accuracy-with-increasing-sequencing-cycles/) for more detailed explaination as to why quality drops with increasing sequencing cycles.

> [A video FastQC walkthrough created by FastQC developers](https://www.youtube.com/watch?v=bz93ReOv87Y "FastQC video") 

## Quality Trimming using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic "Trimmomatic Homepage")
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

Now we will run Trimmomatic on these raw data to remove low quality reads as well as adapters. 

>i. Get an interactive cluster node to start running programs


>ii. Create these output directories in your day1_morn folder to save trimmomatic results

```
mkdir Rush_KPC_266_trimmomatic_results
```

>iii. Load latest version of java and try to invoke trimmomatic from command line.

```
module load lsa java/1.8.0

java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar –h
```

Pending: explaining parameters and its default value. Adapter file. Changing only SLIDINGWINDOW parameter from default 4:15 to 4:20 for raw reads.

>iv. Run the below trimmomatic commands on raw reads.

```
time java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_266_1_combine.fastq.gz Rush_KPC_266_2_combine.fastq.gz Rush_KPC_266_trimmomatic_results/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results/forward_unpaired.fq.gz Rush_KPC_266_trimmomatic_results/reverse_paired.fq.gz Rush_KPC_266_trimmomatic_results/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:40 HEADCROP:0
```

>v. Now create new directories in day1_morn folder and Run FastQC on these trimmomatic results.

```
mkdir Rush_KPC_266_FastQC_results/after_trimmomatic

fastqc -o Rush_KPC_266_FastQC_results/after_trimmomatic/ Rush_KPC_266_trimmomatic_results/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results/reverse_paired.fq.gz --extract
scp username@flux-xfer.engin.umich.edu:/scratch/micro612w16_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/after_trimmomatic/*.html /path-to-local-directory/
```

`screenshots explaination`
`before trimmomatic and after trimmomatic explanation: How quality and overrepresented sequences red cross signal disappeared after running trimmomatic`

If you notice the per base sequence content graph, the head bases(~9 bp) are slightly imbalanced. Each nucleotide content should run parallel to each other. This is not very bad but you can fix this by trimming these imbalanced head bases using HEADCROP:9 flag in the above command.

>vi. Lets Run trimmomatic again with headcrop 9 and save it in a different directory called Rush_KPC_266_trimmomatic_results_with_headcrop/

```
mkdir Rush_KPC_266_trimmomatic_results_with_headcrop/

time java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_266_1_combine.fastq.gz Rush_KPC_266_2_combine.fastq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/forward_unpaired.fq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/reverse_paired.fq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:40 HEADCROP:9
```

>vii. Run FastQC on updated trimmomatic results with headcrop and check report on your local computer

```
mkdir Rush_KPC_266_FastQC_results/after_trimmomatic_headcrop/
fastqc -o Rush_KPC_266_FastQC_results/after_trimmomatic_headcrop/ --extract -f fastq Rush_KPC_266_trimmomatic_results_with_headcrop/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/reverse_paired.fq.gz
scp username@flux-xfer.engin.umich.edu:/scratch/micro612w16_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/after_trimmomatic_headcrop/*.html /path-to-local-directory/
```

Notice the per base sequence content graph in report changed to just warning from red cross sign.
`screenshot explanation`


[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

# Day 1 Afternoon

Read Mapping is one of the most common Bioinformatics operation that needs to be carried out on NGS data. Reads are generally generally mapped to a reference genome sequence or closely related genome. There are number of tools that can map reads to a reference genome and differ from each other in algorithm, speed and accuracy. Most of these tools work first builds an index of the reference sequence which works like a dictionary for fast search/lookup and then calling an alignment algorithm which uses these index to align short read sequences against the reference. These alignment has a vast number of uses ranging from Variant/SNP calling, Coverage estimation and expression analysis.

## Read Mapping
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

**1. Copy day1_after directory from shared data directory in your home directory.**

```
cp -r /scratch/micro612w16_fluxod/shared/data/day1_after ./
cd day1_after/
mkdir Rush_KPC_266_varcall_result
```

We will be using the trimmed clean reads that were obtained after Trimmmatic on raw reads.

**2. Map your reads against a finished reference genome using [BWA](http://bio-bwa.sourceforge.net/bwa.shtml "BWA manual")**

BWA is one of the several and a very good example of read mappers that are based on Burrows-Wheeler transform algorithm. If you feel like challenging yourselves, you can read BWA paper [here](http://bioinformatics.oxfordjournals.org/content/25/14/1754.short)

>i. To create BWA index of Reference fasta file, you need to run the following command.

path to reference fasta: /scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta

> No need to create index, We have already created an index for you using this command.
```
bwa index /scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta
```

>ii. Align reads to reference and output into SAM file

Now lets align both left and right end reads to our reference using BWA alignment algorithm 'mem' which is one of the three algorithms that is fast and works on mate paired end reads. For other algorithms, you can refer to BWA [manual](http://bio-bwa.sourceforge.net/bwa.shtml "BWA manual")

```
bwa mem -M -R "@RG\tID:96\tSM:Rush_KPC_266_1_combine.fastq.gz\tLB:1\tPL:Illumina" -t 8 /scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta forward_paired.fq.gz reverse_paired.fq.gz > Rush_KPC_266__aln.sam
```

> -R readgroup parameter; what does it say? screenshot.

**3. SAM/BAM manipulation and variant calling using [Samtools](http://www.htslib.org/doc/samtools.html "Samtools Manual")**

>i. Change directory to results folder:

```
cd Rush_KPC_266_varcall_result
```

>ii. Convert SAM to BAM using SAMTOOLS:

```
samtools view -Sb Rush_KPC_266__aln.sam > Rush_KPC_266__aln.bam
```

>iii. Sort BAM file using SAMTOOLS:

```
samtools sort Rush_KPC_266__aln.bam Rush_KPC_266__aln_sort
```

**4. Mark duplicates(PCR optical duplicates) and remove them using [PICARD](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates "Picard MarkDuplicates")**

>i. Create a dictionary for reference fasta file required by PICARD(If KPNIH1.dict doesn’t exist. Ignore if already exists.).
 
```
java -jar /scratch/micro612w16_fluxod/shared/bin/picard-tools-1.130/picard.jar CreateSequenceDictionary REFERENCE=/scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta OUTPUT=/scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.dict
```

>ii. Run PICARD for removing duplicates.

```
java -jar /scratch/micro612w16_fluxod/shared/bin/picard-tools-1.130/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=Rush_KPC_266__aln_sort.bam OUTPUT= Rush_KPC_266__aln_marked.bam METRICS_FILE=Rush_KPC_266__markduplicates_metrics
CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
```

>iii. Sort these marked BAM file again (for downstream compatibility)

```
samtools sort Rush_KPC_266__aln_marked.bam Rush_KPC_266__aln_sort
```

>iv. Index these marked bam file using SAMTOOLS

```
samtools index Rush_KPC_266__aln_sort.bam
```

## Variant Calling and Filteration
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

**1. Call variants using [samtools](http://www.htslib.org/doc/samtools.html "samtools manual") mpileup and [bcftools](https://samtools.github.io/bcftools/bcftools.html "bcftools")**

```
samtools mpileup -ug -f /scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta Rush_KPC_266__aln_sort.bam | bcftools call -O v -v -c -o Rush_KPC_266__aln_mpileup_raw.vcf
```

> -g generates genotype likelihood in bcf format   
> -c call samtools consensus caller

**2. Variant filtering and processed file generation using GATK and vcftools**

>i. Variant filtering using [GATK](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php "GATK Variant Filteration"):

Now lets filter these variants based on their quality using GATK.

```
java -jar /scratch/micro612w16_fluxod/shared/bin/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T VariantFiltration -R
/home2/apirani/bin/reference/KPNIH1/KPNIH1.fasta -o Rush_KPC_266__filter_gatk.vcf --variant Rush_KPC_266__aln_mpileup_raw.vcf --filterExpression "FQ < 0.025 && MQ > 50 && QUAL > 100 && DP > 15" --filterName pass_filter
```

Open the file Rush_KPC_266__filter_gatk.vcf and have a look at 7th column. Take a glance at the rows with 'pass_filter' in its 7 column. 

```
grep 'pass_filter' Rush_KPC_266__filter_gatk.vcf | head
```

The quality filters that we generally employ to get high quality variants are:
FQ: Consensus Quality
MQ: Mapping quality.
DP: Depth of reads supporting that variant.
QUAL: Overall quality of the variant called.

But going forward you can use other parameters as well based on your requirements.
More Info on VCF format and parameter specifications can be found [here](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0ahUKEwjzkcSP4MfLAhUDyYMKHU3yDwMQFggjMAA&url=https%3A%2F%2Fsamtools.github.io%2Fhts-specs%2FVCFv4.2.pdf&usg=AFQjCNGFka33WgRmvOfOfp4nSaCzkV95HA&sig2=6Xb3XDaZfghadZfcnnPQxw&cad=rja "VCF format Specs.")

>ii. Remove indels and keep only SNPS that passed our filter criteria using [vcftools](http://vcftools.sourceforge.net/man_latest.html vcftools manual):

In most outbreak studies, we are trying to find how these outbreak samples differ and evolve from the reference/index genome. Many tools that carry out such type of phylogenetic analysis requires a consensus sequences containing only variant calls(SNPs). Though there are few tools that take into consideration both SNPs and Indels and give a greater resolution.
Now we will try to construct a consensus sequence using only SNP calls.

vcftools is a program package that is especially written to work with vcf file formats. It thus saves your precious time by making available all the common operations with a single command.
Lets remove indels from our final vcf file and keep only variants that passed the filters.

```
vcftools --vcf Rush_KPC_266__filter_gatk.vcf --keep-filtered pass_filter --remove-indels --recode --recode-INFO-all --out
Rush_KPC_266__filter_onlysnp 
```

Notice the details that were printed out in STDOUT.

>iii. Generate Consensus fasta file from filtered variants using vcftools:

A consensus fasta sequence will contain alleles from reference sequence at positions where no variants were observed and variants that were observed at positions described in vcf file.

```
bgzip Rush_KPC_266__filter_onlysnp.recode.vcf
tabix Rush_KPC_266__filter_onlysnp.recode.vcf.gz
cat /scratch/micro612w16_fluxod/shared/bin/reference/KPNIH1/KPNIH1.fasta | vcf-consensus Rush_KPC_266__filter_onlysnp.recode.vcf.gz > Rush_KPC_266__consensus.fa
```

Check the fasta header and change it using sed.

```
head -n1 Rush_KPC_266__consensus.fa
sed -i 's/>.*/>Rush_KPC_266_/g' Rush_KPC_266__consensus.fa 
```

**3. Variant Annotation using snpEff**

Variant annotation is one of the crucial steps in any Variant Calling Pipeline. Most of the variant annotation tools creates their own database or an external one to assign function and predicts the effect of variants on genes. We will try to touch base on some basic steps of annotating variants in our vcf file using snpEff. 

>i. Check snpEff internal database for your reference genome:

```     
java -jar /scratch/micro612w16_fluxod/shared/bin/snpEff/snpEff.jar databases | grep 'kpnih1'
```
Note down the genome id for your reference genome KPNIH1. In this case: GCA_000281535.2.29

>ii. Change the Chromosome name in vcf file to ‘Chromosome’ for snpEff reference database compatibility. 

```
sed -i 's/gi.*|/Chromosome/g' Rush_KPC_266__filter_gatk.vcf
```
>iii. Run snpEff for variant annotation.

```
java -jar snpEff.jar -v GCA_000281535.2.29 Rush_KPC_266__filter_gatk.vcf > Rush_KPC_266__filter_gatk_ann.vcf -htmlStats Rush_KPC_266__filter_gatk_stats
```

Check the STDOUT printed which contains some useful details. Explain!

Lets go through the ANN field added after annotation step.

```
grep '^Chromosome' Rush_KPC_266__filter_gatk_ann.vcf | head -n1
```

Explain ANN field!!!

**4. Generate Statistics report using samtools, vcftools and qualimap**

Lets try to get some statistics about various outputs that were created using the above steps and check if everything makes sense.

>i. Reads Alignment statistics:
 
```
samtools flagstat Rush_KPC_266__aln.bam > Rush_KPC_266__alignment_stats
```

ii. VCF statistics:  

```
bgzip Rush_KPC_266__aln_mpileup_raw.vcf   
tabix Rush_KPC_266__aln_mpileup_raw.vcf.gz  
vcf-stats Rush_KPC_266__aln_mpileup_raw.vcf.gz > Rush_KPC_266__raw_vcf_stats
```

iii. Qualimap report of BAM coverage:

``` 
qualimap bamqc -bam Rush_KPC_266__aln_sort.bam -outdir ./ -outfile Rush_KPC_266__report.pdf -outformat pdf 
```

Open the pdf report in your local system.

```
scp username@flux-xfer.engin.umich.edu:/scratch/micro612w16_fluxod/username/day1_after/Rush_KPC_266_varcall_result/Rush_KPC_266__report.pdf /path-to-local-directory/
```

Check the Chromosome stats table.
Check the coverage across reference and Mapping quality across the reference.

## Visualize BAM and VCF files in IGV or ACT(Screenshots Pending)
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

> Required Input files: KPNIH1 reference fasta and genbank file, Rush_KPC_266__aln_marked.bam and Rush_KPC_266__aln_marked.bai, Rush_KPC_266__aln_mpileup_raw.vcf/Rush_KPC_266__filter_onlysnp.recode.vcf/Rush_KPC_266__filter_gatk_ann.vcf

```
screenshots explanation here
```

[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)


# Day 2 Morning
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

On day 1 we worked through a pipeline to map short-read data to a pre-existing assembly and identify single-nucleotide variants (SNVs) and small insertion/deletions. However, what this sort of analysis misses is the existence of sequence that is not present in your reference. Today we will tackle this issue by assembling our short reads into larger sequences, which we will then analyze to characterize the functions unique to our sequenced genome.   

Execute the following command to copy files for this morning’s exercises to your scratch directory: 

```
cd /scratch/micro612w16_fluxod/username
cp -r /scratch/micro612w16_fluxod/shared/data/day2_morn ./
```

## Genome Assembly using Spades Pipeline
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

There is a wide range of tools available for assembly of microbial genomes. These assemblers fall in to two general algorithmic categories, which you can learn more about [here](Ask Evan for the link?). In the end, most assemblers will perform well on microbial genomes, unless there is unusually high GC-content or an over-abundance of repetitive sequences, both of which make accurate assembly difficult. 

Here we will use the Spades assembler with default parameters. Because genome assembly is a computationally intensive process, we will submit our assembly jobs to the cluster, and move ahead with some pre-assembled genomes, while your assemblies are running. 

>i. Create directory to hold your assembly output.

create a new directory for the spades output in your day2_morn folder
```
cd /scratch/micro612w16_fluxod/username/day2_morn
mkdir Rush_KPC_266_assembly_result
```

>ii. Test out Spades to make sure its in your path

To make sure that your paths are set up correctly, try running Spades with the –h (help) flag, which should produce usage instruction.

```
spades.py -h     
OR    
python spades.py –h	
(In case, Invoke doesn’t work, check python path supplied in spades.py script or invoke using python interpreter)
```

>iii. Submit a cluster job to assemble 

Open the spades.PBS with nano and add the following spades command to the bottom of the file. Dont forget to change 'username' with your unique id. Also Change JOBNAME and EMAIL_ADDRESS in pbs script accordingly.

```
python /scratch/micro612w16_fluxod/shared/bin/Spades/bin/spades.py --pe1-1 /scratch/micro612w16_fluxod/apirani/day2_morn/forward_paired.fq.gz --pe1-2 /scratch/micro612w16_fluxod/apirani/day2_morn/reverse_paired.fq.gz --pe1-s /scratch/micro612w16_fluxod/apirani/day2_morn/forward_unpaired.fq.gz --pe1-s /scratch/micro612w16_fluxod/apirani/day2_morn/reverse_unpaired.fq.gz -o /scratch/micro612w16_fluxod/apirani/day2_morn/Rush_KPC_266_assembly_result/ --careful
```

>iv. Submit your job to the cluster with qsub

```
qsub spades.PBS
```

>v. Verify that your job is in the queue with the qstat

```
qstat –u username 
```

## Assembly evaluation using QUAST
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

The output of an assembler is a set of contigs (contiguous sequences), that are composed of the short reads that we fed in. Once we have an assembly we want to evaluate how good it is. This is somewhat qualitative, but there are some standard metrics that people use to quantify the quality of their assembly. Useful metrics include: i) number of contigs (the fewer the better), ii) N50 (the minimum contig size that at least 50% of your assembly belongs, the bigger the better). In general you want your assembly to be less than 200 contigs and have an N50 greater than 50 Kb, although these numbers of highly dependent on the properties of the assembled genome. 

To evaluate some example assemblies we will use the tool quast. Quast produces a series of metrics describing the quality of your genome assemblies. 

>i. Run quast on a set of previously generated assemblies

```
quast.py -o quast sample_264_contigs.fasta sample_266_contigs.fasta
```

>ii. Explore quast output

QUAST creates output in various format. Now lets check the report.txt for assembly statistics. Open report.txt using nano.

```
nano quast/report.txt
```
Check the difference between each assembly statistics. 

## Compare assembly to reference genome and Post-assembly genome improvement
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

Now that we feel confident in our assembly, lets compare it to our reference to see if we can identify any large insertions/deletions. To do this we will use a graphical to call Artemis Comparison Tool (ACT). To do this we need to first align our genome assembly to our reference. We will accomplish this using command-line BLAST.

>i. Align unordered contigs to reference

Create a BLAST database from your reference genome using the makeblastdb command.

```
makeblastdb -in KPNIH1.fasta -dbtype nucl -out KPNIH1.fasta
```

>ii. Stitch together your contigs into a single sequence

```
echo ">sample_266_contigs_concat" > sample_266_contigs_concat.fasta 
grep -v ">" sample_266_contigs.fasta >> sample_266_contigs_concat.fasta 
```

BLAST your stitched together contigs against your reference. The input parameters are: 1) query sequences (-query sample_266_contigs_concat.fasta), 2) the database to search against (-db KPNIH1.fasta), 3) the name of a file to store your results (-out blastn_results), 4) output format (-outfmt 6), 6) e-value cutoff (-evalue 1e-20)

```
blastn -outfmt 6 -evalue 1e-20 -db KPNIH1.fasta -query sample_266_contigs_concat.fasta -out concat_comp.blast
```

>ii. Use ACT to compare stitched together contigs to reference

```
cd /scratch/micro612w16_fluxod/username/day2_morn
mkdir ACT_contig_comparison 
cp KPNIH.gb KPNIH1.fasta concat_comp.blast ACT_contig_comparison/

Use scp to get sequences and BLAST alignments onto your laptop 
scp -r username@flux-xfer.engin.umich.edu:/scratch/micro612w16_fluxod/username/day2_morn/ACT_contig_comparison/ /path-to-local-directory/
```

>iii. Read these Input files in ACT_contig_comparison into ACT

```
Go to File -> open 
Sequence file 1 = KPNIH.gb
Comparison file 1  = concat_comp_blast 
Sequence file 2  = sample_266_contigs_concat.fasta
```

> Notice that it a complete mess!!!! The reason is that the contigs are in random order, so it is very difficult to visually compare to the reference. 

iv. Run abacas to orient contigs to reference

To orient our contigs relative to the reference we will use a tool called abacas. Abacas aligns contigs to a reference genome and then stitches them together to form a “pseudo-chromosome”. Go back to flux and into the directory where the assembly are located.

```
cd /scratch/micro612w16_fluxod/username/day2_morn/
```

Run abacas using the input parameters: 
1) your reference sequence (-r KPNIH.fasta), 
2) your contig file (-q sample_266_contigs.fasta), 
3) the program to use to align contigs to reference (-p nucmer), 
4) append unmapped contigs to end of file (-b), 
5) use default nucmer parameters (-d), 
6) append contigs into pseudo-chromosome (-a), 
7) the prefix for your output files (–o sample_266_contigs_ordered) 

```
perl abacas.1.3.1.pl -r KPNIH1.fasta -q sample_266_contigs.fasta -p nucmer -b -d -a -o sample_266_contigs_ordered
```

v. Use ACT to view contig alignment to reference genome

> Use scp to get ordered fasta sequence and .cruch file onto your laptop 
scp username@flux-xfer.engin.umich.edu:/scratch/micro612w16_fluxod/username/day2_morn/sample_266_contigs_ordered* /path-to-local-ACT_contig_comparison-directory/
```

> Read files into ACT

```
Go to File -> open 
Sequence file 1 = KPNIH.gb 
Comparison file 1  = sample_266_contigs_ordered.crunch 
Sequence file 2  = sample_266_contigs_ordered.fasta
```

> Notice that the alignment is totally beautiful now!!! Scan through the alignment and play with ACT features to look at genes present in reference but not in assembly 

## Map reads to the final ordered assembly(To do or not to do!)
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

## Genome Annotation
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

**Identify protein-coding genes with Prokka**

From our ACT comparison of our assembly and the reference we can clearly see that there is unique sequence in our assembly. However, we still don’t know what that sequence encodes! To try to get some insight into the sorts of genes unique to our assembly we will run a genome annotation pipeline called Prokka. Prokka works by first running denovo gene prediction algorithms to identify protein coding genes and tRNA genes. Next, for protein coding genes Prokka runs a series of comparisons against databases of annotated genes to generate putative annotations for your genome. 

>i. Run Prokka on assembly

Load modules required for Prokka

```
module load med perl-modules prokka 
prokka –setupdb
```

Execute Prokka on your ordered assembly 

```
cd /scratch/micro612w16_fluxod/username/day2_morn/
mkdir sample_266_prokka 
prokka -kingdom Bacteria -outdir /scratch/micro612w16_fluxod/username/day2_morn/sample_266_prokka -force -prefix sample_266 sample_266_contigs_ordered.fasta

Use scp to get Prokka annotated genome on your laptop

scp -r username@flux-xfer.engin.umich.edu:/scratch/micro612w16_fluxod/username/day2_morn/sample_266_prokka/ /path-to-local-ACT_contig_comparison-directory/
```

>ii. Reload comparison into ACT now that we’ve annotated the un-annotated!

Read files into ACT
```
Go to File -> open
Sequence file 1  = KPNIH.gb 
Comparison file 1  = sample_266_contigs_ordered.crunch 
Sequence file 2  = sample_266_contigs_ordered.gbf
```

>Play around with ACT to see what types of genes are unique to sample 266!!! 

# Day 2 Afternoon
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

## High-throughput BLAST and pan-genome analysis

This morning we learned how to perform basic genome annotation and comparison using Prokka and ACT. Now we will up the ante and do some more sophisticated comparative genomics analyses! 
First, we will create custom BLAST databases to identify specific antibiotic resistance genes of interest in a set of genomes. 
Second, we will use the large-scale BLAST-based tool LS-BSR to identify the complete antibiotic resistome in our genomes. 
Third, we will move beyond antibiotic resistance, and look at the complete set of protein coding genes in our input genomes. 
Finally, we will go back to ACT to understand the sorts of genomic rearrangements underlying observed variation in gene content.  

For these exercises we will be looking at four closely related Acinetobacter baumannii strains. However, despite being closely related, these genomes have major differences in gene content, as A. baumannii has a notoriously flexible genome! In fact, in large part due to its genomic flexibility, A. baumannii has transitioned from a harmless environmental contaminant to a pan-resistant super-bug in a matter of a few decades. If you are interested in learning more, check out this nature [review](http://www.nature.com/nrmicro/journal/v5/n12/abs/nrmicro1789.html) or [this](http://www.pnas.org/content/108/33/13758.abstract) paper, I published a few years back analyzing the very same genomes you are working with.

Execute the following command to copy files for this afternoon’s exercises to your scratch directory:

```  
cd /scratch/micro612w16_fluxod/username
cp -r /scratch/micro612w16_fluxod/shared/data/day2_after/ ./
```

##1. Determine which genomes contain beta-lactamase genes
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

Before comparing full genomic content, lets start by looking for the presence of particular genes of interest. A. baumannii harbors an arsenal of resistance genes, and it would be interesting to know how particular resistance families vary among our 4 genomes. To accomplish this we will use the antibiotic resistance database (ARDB). In particular, we are going to extract a set of genes from ARDB that we are interested in probing our genomes for, and create a custom BLAST database to compare against. i. Get beta-lactamase genes from ARDB database

>i. Run the custom perl script filter_fasta_file.pl to extract genes annotated as beta-lactamases from the full ARDB fasta file. 

The script takes as input: 
1) the ARDB database (resisGenes.pfasta), 
2) a file containing terms to search the database for (fasta_file_keys) and 
3) an output file to contain the subset of sequences that match the text your searching for (ardb_beta_lactam_genes.pfasta).

```
cd scratch/micro612w16_fluxod/username/day2_after
perl filter_fasta_file.pl resisGenes.pfasta fasta_file_keys ardb_beta_lactam_genes.pfasta
```

>ii. Build BLAST database from fasta file

Run formatdb on the file of beta-lactamases to create a BLAST database. 
formatdb takes as input: 
1) a fasta file of protein or nucleotide sequences (ardb_beta_lactam_genes.pfasta) and 
2) a flag indicating whether to construct a protein or nucleotide database (in this case protein/ -p T).

```
formatdb -i ardb_beta_lactam_genes.pfasta -p T
```

>iii. BLAST A. baumannii proteins against our custom beta-lactamase database

Run BLAST! 

The input parameters are: 

1) the type of blast to use (-p blastp), 
2) query sequences (-i Abau_all.pfasta), 
3) the database to search against (-d ardb_beta_lactam_genes.pfasta), 
4) the name of a file to store your results (-o bl_blastp_results), 
5) output format (-m 8), 
6) e-value cutoff (-e 1e-20), 
7) number of database sequences to return (-v 1) and 
8) number of database sequences to show alignment for (-b 1).

```
blastall -p blastp -i Abau_all.pfasta -d ardb_beta_lactam_genes.pfasta -o bl_blastp_results -m 8 -e 1e-20 -v 1 -b 1
```

Use less to look at bl_blastp_results.
Experiment with the –m parameter, which controls different output formats that BLAST can produce. 


>iv. Repeat steps i-iii for a different resistance gene class

Use nano to change fasta_file_keys to contain phrase you’d like to search for (e.g. acetyltransferase, carbapenemase)

Run filter_fasta_file.pl to extract genes annotated with your resistance of interest (ROI) from the full ARDB fasta file

```
perl filter_fasta_file.pl resisGenes.pfasta fasta_file_keys ardb_ROI_genes.pfasta
```

Create new BLAST database with formatdb

```
formatdb -i ardb_ROI_genes.pfasta -p T
```

BLAST!

```
blastall -p blastp -i Abau_all.pfasta -d ardb_ROI_genes.pfasta -o bl_blastp_results -m 8 -e 1e-20 -v 1 -b 1
```

##2. Identification of antibiotic resistance genes with LS-BSR and the ARDB database
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

Next, instead of looking at resistance classes one at a time, lets look at them all in one shot! To do this we will use LS-BSR, which essentially is just a wrapper for doing the same sort of BLASTing we just did in the previous step. BSR stands for BLAST Score Ratio, which refers to what the output is. In particular, for each query gene LS-BSR returns the ratio between: 1) the BLAST score of best hit in target genome and 2) BLAST score of query gene against itself. So, the output is a query x target genome matrix, where the values are between 0 and 1, and indicate the strength of a given queries BLAST hit in the target genome. 

>i. Load modules required for LS-BSR

```
module load med 
module load sph 
module load lsa 
module load usearch 
module load python/2.7.3 
module load biopython 
module load ls-bsr 
module load prodigal
```

>ii. Create a non-redundant list of resistance genes

There is a lot of redundancy in the ARDB (e.g. lots of closely related genes), which would make the output difficult to sort through. Here, we use usearch to select representatives from the database and create a non-redundant gene set! 

We are running usearch with the following parameters: 
1) the clustering algorithm (-cluster_fast), 
2) the files of sequences you want to cluster (resisGenes.pep), 
3) the minimum sequence identity to be included in an existing cluster (-id 0.8), 
4) an output fasta file with reperesentatives (centroids) of each sequence cluster (-centroids resisGenes_nr.pep) and 
5) an output file describing the results of the clustering (-uc resisGenes.uc).

```
cd scratch/micro612w16_fluxod/username/day2_after
usearch -cluster_fast resisGenes.pep -id 0.8 -centroids resisGenes_nr.pep -uc resisGenes.uc
```

>iii. Run LS-BSR

LS-BSR is pretty intensive, so we want to get an interactive node to run this

```
qsub -I -V -l nodes=1:ppn=1,mem=4000mb,walltime=7:00:00:00 -q fluxod -l qos=flux -A micro612w16_fluxod
```

Run LS-BSR (it will take a few minutes)! 

The input parameters are: a directory with your genomes (-d Abau_genomes) and a fasta file of query genes (-g resisGenes_nr.pep)

```
cd scratch/micro612w16_fluxod/username/day2_after 
python /home/software/rhel6/med/python-libs/ls-bsr/1.0/LS-BSR-master/ls_bsr.py -d Abau_genomes/ -g resisGenes_nr.pep
```

>iv. Download LS-BSR output matrix to your own computer for analysis in R

Use sftp to get LS-BSR output onto your laptop

```
cd ~/Desktop (or wherever your desktop is) 
mkdir LS-BSR_resistome 
cd LS-BSR_resistome 
sftp –r username@flux-login.engin.umich.edu cd /scratch/micro612w16_fluxod/username/day2_after 
get bsr_matrix_values.txt
```

Fire up RStudio and read the matrix:

```
bsr_mat = read.table('bsr_matrix_values.txt', sep = "\t", row.names = 1, header = TRUE, quote = "")
```

Use head, str, dim, etc. to explore the matrix you read in

v. Make a heatmap of all the LS-BSR results

Install and load the R library "heatmap3"

Make a heatmap of the complete LS-BSR matrix. Check out the help file to see what the input parameters do, and behold the plethora of other options to customize your heatmaps!

```
heatmap3(bsr_mat, , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5)
```


>vi. Subset LS-BSR data to only include genes present in at least one genome 

From the previous step you should have discerned that full LS-BSR matrix is too large to get a useful visualization, so we need to subset it. 
Lets first subset the matrix to focus only on genes present in at least one of our genomes. 
Values in the LS-BSR matrix are between 0 and 1, and represent the sequence identity to the query gene. 
We will arbitrarily say that if a protein have a BLAST score ratio of less then 0.5, then its absent.

```
bsr_mat_subset = bsr_mat[rowSums(bsr_mat > 0.5) > 0,]
```

Make a heatmap of your subset (much better!)

```
heatmap3(bsr_mat_subset, , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5)
```

>vii. Determine the total number of resistance genes present in each genome

We use colSums to count the number of genes with greater than 50% identity to the query

```
colSums(bsr_mat > 0.5)
```

How does the total number of genes vary by altering the percent identity threshold?

>viii. Determine the total number of bla genes in each genome

Next, we will use grepl to pull out genes of interest

```
bla_bsr_mat = bsr_mat[grepl('beta-lactamase', row.names(bsr_mat)) ,]
```

Print out to screen and make a heatmap to explore

>ix. Subset the full matrix to look at genes that are present in only one genome

Get genes present in only one genome

```
unique_bsr_mat = bsr_mat[rowSums(bsr_mat > .5) == 1,]
```

Print out to screen and make a heatmap to explore

##3. Perform pan-genome analysis with LS-BSR
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

As a final BLASTing exercise we will use LS-BSR to explore the pan-genome of our A. baumannii. The pan-genome is just a fancy term for the full complement of genes in a set of genomes. 
The way LS-BSR does this is by: 
1) applying prodigal to identify protein coding genes in input genomes, 
2) applying usearch to create non-redundant set of genes and 
3) BLASTing the set of non-redundant genes against the genomes.

>i. Get pan-genome matrix and transfer annotation

Make sure you are on an interactive node, as this will be even more computationally intensive!

```
qsub -I -V -l nodes=1:ppn=1,mem=4000mb,walltime=7:00:00:00 -q fluxod -l qos=flux -A micro612w16_fluxod
```

Run LS-BSR! The –u parameter is just a path to where usearch lives on flux.

```
cd scratch/micro612w16_fluxod/username/day2_after

python /home/software/rhel6/med/python-libs/ls-bsr/1.0/LS-BSR-master/ls_bsr.py -d ../Abau_genomes/ -u /home/software/rhel6/sph/usearch/7.0.1001/bin/usearch7.0.1001_i86linux32
```

Run the custom perl script transfer_annotations.pl to add annotations to your BSR matrix. The output of this script will be bsr_matrix_values_annot.txt

```
perl transfer_annotations.pl Abau_ECII_PC.fasta Abau_ECII_PC.NR.annot bsr_matrix_values.txt consensus.fasta
```

>ii. Read matrix into R and create heatmap

Use sftp to get LS-BSR output onto your laptop

```
cd ~/Desktop (or wherever your desktop is) 
mkdir LS-BSR_pan_genome 
cd LS-BSR_pan_genome 
sftp –r username@flux-login.engin.umich.edu cd /scratch/micro612w16_fluxod/username/day2_after get bsr_matrix_values_annot.txt
```

Fire up RStudio and read the matrix in

```
bsr_mat_PG = read.table('bsr_matrix_values_annot.txt', sep = "\t", row.names = 1, header = TRUE, quote = "")
```

Use head, str, dim, etc. to explore the matrix you read in
Make a heatmap for the full matrix

```
heatmap3(as.matrix(bsr_mat_PG), , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5)
```

Make a heatmap for variable genes (present in at least one, but not all of the genomes

```
bsr_mat_PG_subset = bsr_mat_PG[rowSums(bsr_mat_PG > 0.4) > 0 & rowSums(bsr_mat_PG > 0.4) < 4 ,] heatmap3(as.matrix(bsr_mat_PG_subset), , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5)
```

>iii. Which genomes are most closely related based upon shared gene content?

We will use the outer function to determine the number of genes shared by each pair of genomes. 
Here we are arbitrarily deciding that a gene is present if the BSR is greater than 0.4. 
Look at the help page for outer to gain additional insight into how this is working.

```
outer(1:4,1:4, FUN = Vectorize(function(x,y){sum(bsr_mat_PG_subset[,x] > 0.4 & bsr_mat_PG_subset[,y] > 0.4)}))
```

>iv. What is the size of the core genome?

Lets first get an overview of how many genes are present in different numbers of genomes (0, 1, 2, 3 or 4) by plotting a histogram. Here, we combine hist with rowSums to accomplish this.

```
hist(rowSums(bsr_mat_PG > 0.4))
```

Next, lets figure out how big the core genome is (e.g. how many genes are common to all of our genomes)?

```
sum(rowSums(bsr_mat_PG > 0.4) == 4)
```

>v. What is the size of the accessory genome?

Lets use a similar approach to determine the size of the accessory genome (e.g. those genes present in only a subset of our genomes).

```
sum(rowSums(bsr_mat_PG > 0.4) < 4 & rowSums(bsr_mat_PG > 0.4) > 0)
```

>vi. What types of genes are unique to a given genome?

So far we have quantified the core and accessory genome, now lets see if we can get an idea of what types of genes are core vs. accessory. Lets start by looking at those genes present in only a single genome. What do you notice about these genes?

```
row.names(bsr_mat_PG[rowSums(bsr_mat_PG > 0.4) == 1,])
```

vii. What is the number of hypothetical genes in core vs. accessory genome?

Looking at unqiue genes we see that many are annotated as “hypothetical”, indicating that the sequence looks like a gene, but has no detectable homology with a functionally characterized gene. Determine the fraction of “hypothetical” genes in unique vs. core. Why does this make sense?

```
sum(grepl("hypothetical" , row.names(bsr_mat_PG[rowSums(bsr_mat_PG > 0.4) == 1,]))) / sum(rowSums(bsr_mat_PG > 0.4) == 1)

sum(grepl("hypothetical" , row.names(bsr_mat_PG[rowSums(bsr_mat_PG > 0.4) == 4,]))) / sum(rowSums(bsr_mat_PG > 0.4) == 4)
```

**4. Perform genome comparisons with ACT**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

In the previous exercises we were focusing on gene content, but losing the context of the structural variation underlying gene content variation (e.g. large insertions and deletions). 
Here we will use ACT to compare two of our genomes (note that you can use ACT to compare more than two genomes if desired). 

i. Create ACT alignment file with BLAST

As we saw this morning, to compare genomes in ACT we need to use BLAST to create the alignments. We will do this on flux.

```
cd scratch/micro612w16_fluxod/username/day2_after
blastall -p blastn -i ../Abau_genomes/AbauA_genome.fasta -d ../Abau_BLAST_DB/ACICU_genome.fasta -m 8 -e 1e-20 -o AbauA_vs_ACICU.blast
```

>ii. Read in genomes, alignments and annotation files

Use sftp to get ACT files onto your laptop

```
cd ~/Desktop (or wherever your desktop is)
mkdir Abau_ACT cd Abau_ACT 
sftp –r username@flux-login.engin.umich.edu 
cd /scratch/micro612w16_fluxod/username/day2_after get Abau_genomes/AbauA_genome.fasta get Abau_genomes/ACICU_genome.fasta 
get AbauA_vs_ACICU.blast 
get Abau_ACT_files/AbauA_genome_gene.gff 
get Abau_ACT_files/ACICU_genome_gene.gff
```

>iii. Explore genome comparison and features of ACT

Read in genomes and alignment into ACT

```
Go to File -> open 
Sequence file 1  = ACICU_genome.fasta 
Comparison file 1  = AbauA_vs_ACICU.blast
Sequence file 2  = AbauA_genome.fasta
```

Before we use annotation present in genbank files. Here we will use ACT specific annotation files so we get some prettier display (resistance genes = red, transposable elements = bright green)  

```
Go to File -> ACICU_genome.fasta -> Read an entry file = ACICU_genome_gene.gff

Go to File -> AbauA_genome.fasta -> Read an entry file = AbauA_genome_gene.gff
```

Play around in ACT to gain some insight into the sorts of genes present in large insertion/deletion regions. 
See if you can find: 
1) differences in phage content, 
2) membrane biosynthetic gene cluster variation and 
3) antibiotic resistance island variation.


# Day 3 Morning
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

On day 1, we ran through a pipeline to map reads against a reference genome and call variants, but didn’t do much with the variants we identified. Among the most common analyses to perform on a set of variants is to construct phylogenetic trees. Here we will explore different tools for generating and visualizing phylogenetic trees, and also see how recombination can distort phylogenetic signal.

For the first several exercises, we will use the A. baumannii genomes that we worked with yesterday afternoon. 
The backstory on these genomes is that Abau_A, Abau_B and Abau_C are representatives of three clones (as defined by pulsed-field gel electrophoresis - a low-resolution typing method) that were circulating in our hospital. 

One of the goals of our published study was to understand the relationship among these clones to discern whether: 
1) the three clones represent three independent introductions into the hospital or 
2) the three clones originated from a single introduction into the hospital, with subsequent genomic rearrangement leading to the appearance of unique clones. 

The types of phylogenetic analyses you will be performing here are the same types that we used to decipher this mystery.
The other two genomes you will be using are ACICU and AB0057. ACICU is an isolate from a hospital in France, and its close relationship to our isolates makes it a good reference for comparison. AB0057 is a more distantly related isolate that we will utilize as an out-group in our phylogenetic analysis. The utility of an out-group is to help us root our phylogenetic tree, and gain a more nuanced understanding of the relationship among strains.

Execute the following command to copy files for this afternoon’s exercises to your scratch directory:

```
cd /scratch/micro612w16_fluxod/username
cp -r /scratch/micro612w16_fluxod/shared/data/day3_morn ./
```

**1. Perform whole genome alignment with Mauve and convert alignment to other useful formats**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

An alternative approach for identification of variants among genomes is to perform whole genome alignments of assemblies. If the original short read data is unavailable, this might be the only approach available to you. Typically, these programs don’t scale well to large numbers of genomes (e.g. > 100), but they are worth being familiar with. We will use the tool mauve for constructing whole genome alignments of our five A. baumannii genomes.

>i. Perform mauve alignment and transfer xmfa back to flux

Use sftp to get genomes onto your laptop

```
cd ~/Desktop (or wherever your desktop is) 
mkdir Abau_mauve
cd Abau_mauve 
sftp –r username@flux-login.engin.umich.edu 
cd /scratch/micro612w16_fluxod/username/day3_morn 
get Abau_genomes
```

Run mauve to create multiple alignment

i. Open mauve 
ii. File -> align with progressiveMauve 
iii. Click on “Add Sequnce” and add each of the 5 genomes you just downloaded 
iv. Name the output file “mauve_ECII_outgroup” and make sure it is in the directory you created for this exercise 
v. Click Align! 
vi. Wait for Mauve to finish and explore the graphical interface

Use sftp to transfer your alignment back to flux for some processing

```
cd ~/Desktop/Abau_mauve
sftp –r username@flux-login.engin.umich.edu 
cd /scratch/micro612w16_fluxod/username/day3_morn 
put mauve_ECII_outgroup
```
 
>ii. Convert alignment to fasta format

Mauve produces alignments in .xmfa format (use less to see what this looks like), which is not compatible with other programs we want to use. We will use the custom script convert_msa_format.pl to change the alignment to fasta format

```
perl convert_msa_format.pl -i mauve_ECII_outgroup -o mauve_ECII_outgroup.fasta -f fasta -c
```

**2. Perform some DNA sequence comparisons and phylogenetic analysis in ape**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

There are lots of options for phylogenetic analysis. Here, we will use the ape package in R to look at our multiple alignments and construct a tree using the Neighbor Joining method. 

Note that ape has a ton of useful functions for more sophisticated phylogenetic analyses!

>i. Get fasta alignment to your own computer

```
cd ~/Desktop/Abau_mauve
sftp –r username@flux-login.engin.umich.edu 
cd /scratch/micro612w16_fluxod/username/day3_morn 
get mauve_ECII_outgroup.fasta
```

ii. Read alignment into R

Fire up RStudio and install/load ape

Use the read.dna function in ape to read in you multiple alignments. 
Print out the variable to get a summary.

```
abau_msa = read.dna('mauve_ECII_outgroup.fasta', format = "fasta") 
```

>iii. Get variable positions

The DNA object created by read.dna can also be addressed as a matrix, where the columns are positions in the alignment and rows are your sequences. We will next treat our alignment as a matrix, and use apply and colSums to get positions in the alignment that vary among our sequences. Examine these commands in detail to understand how they are working together to give you a logical vector indicating which positions vary in your alignment.

```
abau_msa_bin = apply(abau_msa, 2, FUN = function(x){x == x[1]}) abau_var_pos = colSums(abau_msa_bin) < 5
```

>iv. Get non-gap positions

For our phylogenetic analysis we want to focus on the core genome, so we will next identify positions in the alignment where all our genomes have sequence.

```
non_gap_pos = colSums(as.character(abau_msa) == '-') == 0
```

>v. Count number of variants between sequences

Now that we know which positions in the alignment are core and variable, we can extract these positions and count how many variants there are among our genomes. Do count pairwise variants we will use the dist.dna function in ape. The model parameter indicates that we want to compare sequences by counting differences. Print out the resulting matrix to see how different our genomes are.

```
abau_msa_var = abau_msa[,var_pos & non_gap_pos ]
var_count_matrix = dist.dna(abau_msa_var, model = "N")
```

>vi. Construct phylogenetic tree

Now we are ready to construct our first phylogenetic tree! 

We are going to use the Neighbor Joining algorithm, which takes a matrix of pairwise distances among the input sequences and produces the tree with the minimal total distance. In essence, you can think of this as a distance-based maximum parsimony algorithm, with the advantage being that it runs way faster than if you were to apply a standard maximum parsimony phylogenetic reconstruction.

As a first step we are going to build a more accurate distance matrix, where instead of counting variants, we will measure nucleotide distance using the Jukes-Cantor model of sequence evolution. This is the simplest model of sequence evolution, with a single mutation rate assumed for all types of nucleotide changes.

```
dna_dist_JC = dist.dna(abau_msa, model = "JC")
```

Next, we will use the ape function nj to build our tree from the distance matrix

```
abau_nj_tree = nj(dna_dist_JC)
```

Finally, plot your tree to see how the genomes group.
```
plot(abau_nj_tree)
```

**3. Perform SNP density analysis to discern evidence of recombination**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

An often-overlooked aspect of a proper phylogenetic analysis is to exclude recombinant sequences. Homologous recombination in bacterial genomes is a mode of horizontal transfer, wherein genomic DNA is taken up and swapped in for a homologous sequence. The reason it is critical to account for these recombinant regions is that these horizontally acquired sequences do not represent the phylogenetic history of the strain of interest, but rather in contains information regarding the strain in which the sequence was acquired from. One simple approach for detecting the presence of recombination is to look at the density of variants across a genome. The existence of unusually high or low densities of variants is suggestive that these regions of aberrant density were horizontally acquired. Here we will look at our closely related A. baumannii genomes to see if there is evidence of aberrant variant densities.

>i. Subset sequences to exclude the out-group

For this analysis we want to exclude the out-group, because we are interested in determining whether recombination would hamper our ability to reconstruct the phylogenetic relationship among our closely related set of genomes.  

>Note that the names of the sequences might be different for you, so check that if the command doesn’t work.

```
abau_msa_no_outgroup = abau_msa[c('ACICU_genome.fa/1-3996847','AbauA_genome.fa/1-3953855','AbauC_genome.fa/1-4200364','AbauB_genome.fa/1-4014916'),]
```

>ii. Get variable positions

Next, we will get the variable positions, as before

```
abau_msa_no_outgroup_bin = apply(abau_msa_no_outgroup, 2, FUN = function(x){x == x[1]}) abau_no_outgroup_var_pos = colSums(abau_msa_no_outgroup_bin) < 4
```

>iii. Get non-gap positions

Next, we will get the core positions, as before

```
abau_no_outgroup_non_gap_pos = colSums(as.character(abau_msa_no_outgroup) == '-') == 0
```

>iv. Create overall histogram of SNP density

Finally, create a histogram of SNP density across the genome. Does the density look even, or do you think there might be just a touch of recombination?

```
hist(which(abau_no_outgroup_var_pos & abau_no_outgroup_non_gap_pos), 10000)
```

**4. Perform recombination filtering with gubbins**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

Now that we know there is recombination, we know that we need to filter out the recombinant regions to discern the true phylogenetic relationship among our strains. In fact, this is such an extreme case (~99% of variants of recombinant), that we could be totally misled without filtering recombinant regions. To accomplish this we will use the tool gubbins, which essentially relies on elevated regions of variant density to perform recombination filtering.

>i. Run gubbins on your fasta alignment

Go back on flux and load modules required by gubbins

```
module load python/2.7.3 biopython dendropy reportlab fasttree RAxML fastml/gub gubbins
```

Run gubbins on your fasta formatted alignment

```
cd /scratch/micro612w16_fluxod/username/day3_morn
run_gubbins.py -v -f 50 -o Abau_AB0057 mauve_ECII_outgroup.fasta
```

>ii. Create gubbins output figure

Gubbins produces a series of output files, some of which can be run through another program to produce a visual display of filtered recombinant regions. Run the gubbins_drawer.py script to create a pdf visualization of recombinant regions. 
The inputs are: 
1) the recombination filtered tree created by gubbins (mauve_ECII_outgroup.final_tree.tre),
2) the pdf file to create (mauve_ECII_outgroup.recombination.pdf) and 
3) a .embl representation of recombinant regions (mauve_ECII_outgroup.recombination_predictions.embl).

```
gubbins_drawer.py -t mauve_ECII_outgroup.final_tree.tre -o mauve_ECII_outgroup.recombination.pdf mauve_ECII_outgroup.recombination_predictions.embl
```
>iii. Download and view gubbins figure and filtered tree

Use sftp to get gubbins output files

```
cd ~/Desktop/Abau_mauve
sftp –r username@flux-login.engin.umich.edu 
cd /scratch/micro612w16_fluxod/username/day3_morn 
get mauve_ECII_outgroup.recombination.pdf 
get mauve_ECII_outgroup.final_tree.tre
```

Open up the pdf and observe the recombinant regions filtered out by gubbins. Does it roughly match your expectations based upon your SNP density plots?

Finally, lets look at the recombination-filtered tree to see if this alters our conclusions. 
To view the tree we will use Seaview, which is a multi-purpose tool for: 
1) visualization/construction of multiple alignments and 
2) phylogenetic tree construction. 
Here, we will just use Seaview to view our gubbins tree.

```
In seaview: 
Go to Trees -> import tree (mauve_ECII_outgroup.final_tree.tre) 
To view sub-tree of interest click on “sub-tree” and select the sub-tree excluding the out-group
```


How does the structure look different than the unfiltered tree?

>Note that turning back to the backstory of these isolates, Abau_B and Abau_C were both isolated first from the same patient. So this analysis supports that patient having imported both strains, which likely diverged at a prior hospital at which they resided.

**5. Create annotated publication quality trees with iTOL**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

For the final exercise we will use a different dataset, composed of USA300 methicillin-resistant Staphylococcus aureus genomes. USA300 is a strain of growing concern, as it has been observed to cause infections in both hospitals and in otherwise healthy individuals in the community. An open question is whether there are sub-clades of USA300 in the hospital and the community, or if they are all the same. Here you will create an annotated phylogenetic tree of strains from the community and the hospital, to discern if these form distinct clusters.

>i. Download MRSA genome alignment from flux

Use sftp to get genomes onto your laptop

```
cd ~/Desktop (or wherever your desktop is) 
mkdir MRSA_genomes 
cd MRSA_genomes
sftp –r username@flux-login.engin.umich.edu 
cd /scratch/micro612w16_fluxod/username/day3_morn 
get 2016-3-9_KP_BSI_USA300.fa 
get 2016-3-9_KP_BSI_USA300_iTOL_HA_vs_CA.txt
```

>ii. Look at SNP density for MRSA alignment in R

Before we embark on our phylogenetic analysis, lets look at the SNP density to verify that there is no recombination

```
mrsa_msa = read.dna('2016-3-9_KP_BSI_USA300.fa', format = 'fasta') 
mrsa_msa_bin = apply(mrsa_msa, 2, FUN = function(x){x == x[1]}) 
mrsa_var_pos = colSums(mrsa_msa_bin) < nrow(mrsa_msa_bin) 
hist(which(mrsa_var_pos), 10000)
```

Does it look like there is evidence of recombination?

>iii. Create fasta alignment with only variable positions

Next, lets create a new fasta alignment file containing only the variant positions, as this will be easier to deal with in Seaview

```
write.dna(mrsa_msa[, mrsa_var_pos], file = '2016-3-9_KP_BSI_USA300_var_pos.fa', format = 'fasta')
```

>iv. Read alignment into Seaview and construct Maximum Likelihood tree

In the previous exercise, we used Seaview to look at a pre-existing tree, here we will use Seaview to create a tree from a
multiple sequence alignment 

Read in multiple alignment of variable positions

```
Go to File -> open ('2016-3-9_KP_BSI_USA300_var_pos.fa)
```

Construct maximum likelihood phylogenetic tree with PhyML and default parameters (note, this will take a few minutes)

```
Go to Trees -> PhyML (Select Bootstrap with 20 replicates)
```

Save your tree

```
File -> Save rooted tree
```

Note that in your research it is not a good idea to use these phylogenetic tools completely blind and I strongly encourage embarking on deeper learning yourself, or consulting with an expert before doing an analysis for a publication

v. Read tree into iTOL

```
To make a prettier tree and add annotations we will use iTOL (http://itol.embl.de/). 
Go to http://itol.embl.de/
To load your tree, click on upload, and select the rooted tree you just created in Seaview
```

Explore different visualization options for your tree (e.g. make it circular, show bootstrap values, try collapsing nodes/branches)

Note that you can always reset your tree if you are unhappy with the changes you’ve made

>vi. Add annotations to tree

One of the most powerful features of iTOL is its ability to overlay diverse types of descriptive meta-data on your tree (http://itol.embl.de/help.cgi#datasets). Here, we will overlay our data on whether an isolate was from a community or hospital infection. To do this simply drag-and-drop the annotation file (2016-3-9_KP_BSI_USA300_iTOL_HA_vs_CA.txt) on your tree and voila! 

Do community and hospital isolates cluster together, or are they inter-mixed?




# Day 3 Afternoon
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

## Klebsiella pneumoniae comparative genomic analysis 

To finish up the workshop we are going to go through the process of working up a complete dataset, from start to finish.  This set of genomes originated from a regional outbreak of bla-KPC carrying Klebsiella pneumoniae – one of the most concerning healthcare associated pathogens. 
The goal is to follow up on a previously published epidemiologic analysis (http://cid.oxfordjournals.org/content/53/6/532.abstract), and see if genomics supports prior epidemiologic conclusions and can provide additional insights. We have our genomes, and we know in which regional facility each isolate originated. 

The goal of this exercise is to:
1) process our genomes (QC, variant calling), 
2) perform a phylogenetic analysis and 
3) overlay our meta-data. 

To make this more difficult, the instructions will be much more vague than in previous sessions, and you will be challenged to use what you have learned, both in the past three days and in the prior workshop, to complete this analysis. Hopefully we’ve prepared you to take on the challenge, but remember this is an open book test! 

Feel free to lean on materials from the workshops, manuals of tools and Google (and of course instructors and neighbors). 

Execute the following command to copy files for this afternoon’s exercises to your scratch directory:

```
cd /scratch/micro612w16_fluxod/username
cp –r  /scratch/micro612w16_fluxod/shared/data/day3_after .
```

**1. Perform QC on fastq files**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

On the first morning you ran FastQC to evaluate the quality of a single genome. However, a typical project will include many genomes and you will want to check the quality of all of your samples. From the bash workshop, I hope you can appreciate that you do not want to process 100 genomes by typing 100 commands – rather you want to write a short shell script to do the work for you!


>i. Write a shell script to run FastQC on all fastq files

The fastq files are located in:

```
/scratch/micro612w16_fluxod/shared/data/day3_after_fastq/
```

Rather than copying these to your directory, analyze the files in that directory directly, so everyone doesn’t have to copy 25G to their home directories. 

**HINTS** 
- Your shell script will include a for loop that loops over all of the genomes in the target directory
- The tricky part of this exercise is that each fastq command contains two files (forward and reverse reads). So, you need to take advantage of the fact that the forward and reverse read files both have the same prefix, and you can loop over these prefixes. 
- You should be able to get prefixes by piping the following unix commands: ls, cut, sort, uniq

>ii. Examine output of FastQC to verify that all samples are OK

**2. Examine results of SPANDx pipeline**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

On the afternoon of day 1 we saw how many steps are involved in calling variants relative to a reference genome. However, the same steps are applied to every sample, which makes this very pipeline friendly!  So, you could write your own shell script to string together these commands, or take advantage of one of several published pipelines. Here, we will use the output of the SPANDx pipeline, which takes as input a directory of fastq files and produces core variant and indel calls.

Because it takes a while to run, we have pre-run it for you. Your task will be to sort through the outputs of SPANDx.

>i. Look at overall statistics for variant calling in excel

SPANDx produces an overall summary file of its run that includes:
1) numbers of SNPs/indels, 
2) numbers of filtered SNPs/indels and 
3) average coverage across the reference genome. 

This summary file is in:  Outputs/Single_sample_summary.txt

Use less to look at this file and then apply unix commands to extract and sort individual columns 

**HINTS**
The following unix commands can be used to get sorted lists of coverage and numbers of SNPs/indels: tail, cut, sort

>ii. Look at filtered variants produced by SPANDx in excel

SPANDx also produces a summary file of the variants/indels it identified in the core genome. 

This summary file is: Outputs/Comparative/All_SNPs_annotated.txt 

Use sftp to download this file and view in excel

- View SPANDx manual for interpretation of different columns which can be found [here](https://github.com/dsarov/SPANDx/blob/master/SPANDx%20Manual _v3.1.pdf)
- Back on Flux, use grep to pull SNPs that have HIGH impact
- What types of mutations are predicted to have “HIGH” impact?
- How many genomes do these HIGH impact mutations tend to be present in? How do you interpret this?

**3. Recombination detection and tree generation**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

>i. Plot the distribution of variants across the genome in R

The positions of variants are embedded in the first column of Outputs/Comparative/All_SNPs_annotated.txt, but you have to do some work to isolate them! 

**HINTS**  
- You will need to pipe together two “cut” commands: the first command will use tab as a delimiter and the second will use _. 
- Note that for cut you can specify tab as the delimiter as follows: cut –d$’\t’ and _ as: cut -d ‘_’
- You should redirect the output of your cut commands (a list of SNP positions) to a file called ‘snp_positions.txt’.
- Finally, download this file, read it into R using ‘read.table’ and use ‘hist’ to plot a histogram of the positions
- Do you observe clustering of variants that would be indicative of recombination?

>ii.  Create fasta file of variants from nexus file

SPANDx creates a file of core SNPs in a slightly odd format (transposed nexus). 
This file is called: Outputs/Comparative/Ortho_SNP_matrix.nex

For convenience, apply the custom perl script located in the same directory to convert it to fasta format

```
perl transpose_nex_to_fasta.pl Ortho_SNP_matrix.nex
```

This file Outputs/Comparative/Ortho_SNP_matrix.fasta should now exist

>iii. Create maximum likelihood tree in Seaview

```
Download Ortho_SNP_matrix.fasta to your home computer
Import the file into Seaview and construct a tree using PhyML (100 bootstraps)
Save tree for later analysis
```

**4. Phylogenetic tree annotation and visualization**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

>i. Load the maximum likelihood tree into iTOL

Note that because the out-group is so distantly related it is difficult to make out the structure of the rest of the tree. To remedy this: 
- Click on the KPNIH1 leaf, go to the “tree structure” menu and “delete leaf” 
- Click on the extended branch leading to where KPNIH1 was, go to the “tree structure” menu and click “collapse branch”

>ii. Load the annotation file ‘Rush_KPC_facility_codes_iTOL.txt’ to view the facility of isolation  Play with tree visualization properties to understand how isolates group by facility o Circular vs. normal tree layout o Bootstrap values o Ignoring branch lengths

```
Which facilities appear to have a lot of intra-facility transmission based on grouping of isolates from the same facility? 
Which patient’s infections might have originated from the blue facility?
```

**5. Assessment of genomic deletions**
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

>i. Download genome coverage bed file and load into R

This file is located in: Outputs/Comparative/Bedcov_merge.txt
This file contains information regarding locations in the reference genome that each sequenced genome does and does not map to.
The first 3 columns of the file are:
1) the name of the reference, 
2) the start coordinate of the window and 
3) the end coordinate of the window

The remaining columns are your analyzed genomes, with the values indicating the fraction of the window covered by reads in that genome.

In essence, this file contains information on parts of the reference genome that might have been deleted in one of our sequenced genomes.

After you download this file, read it into R

**HINTS**
- Use the read.table function with the relevant parameters being: header and sep

>ii. Plot heatmap of genome coverage bed file

**HINTS**
- The first 3 columns of the bed file specify the name of the chromosome and the genome coordinates – therefore you want to subset your matrix to not include these columns 
- Use the heatmap3 function to make your heatmap with the following parameters: scale = “none” (keeps original values), Rowv = NA (suppress clustering by rows – why might we not want to cluster by rows for this analysis?)

Note a large genomic deletion among a subset of isolates. Does this deletion fit with the phylogeny from above?

iii. Explore genomic deletion in more detail with ACT

- Use abacus to orient contigs from Rush_KPC_298 to KPNIH 
- Load KPNIH.gb, Rush_KPC_298_ordered and the .crunch alignment into ACT

```
What genes appear to have been lost?
```

# Helpful resources for microbial genomics

1. General Bioinformatics resources
	- [Omictools](http://omictools.com/)
	- [Bioinformatics One-liners by Stephen Turner](https://github.com/stephenturner/oneliners)
	
	
2. Short read processing 
	- [bwa]()
	- [bowtie]()
	- [samtools]()
	- [vcftools]()
	- [bcftools]()
	- [gatk]() 
	- [picard]() 
	- [SPANDx]() 
	- [Snippy]() 
	
3. Genome assembly 
	- [Spades]()
	- [Velvet]() 
	- [Mira]() 
	- [A5]()
	
4. Genome alignment 
	- [Mauve]()
	- [MUMmer]()
	- [Mugsy]() 
	
5. Visualization of genomic data
	- [Artemis]() 
	- [Artemis Comparison Tool]()
	- [IGV]()
	
6. Genome annotation 
	- [Prokka]() 
	- [Blastall]()
	- [LS-BSR]() 
	
7. Phylogenetic tools and resources
	- Visualization 
		- [Seaview]() 
		- [iTOL]() 
		- [Figtree]() 
	- Phylogenetic software 
		- [PAUP]()
		- [RaxML]() 
		- [PhyML]()
		- [BEAST]()
		- [PHYLIP]() 
	- Recombination detection 
		- [Gubbins]()
		- [ClonalFrame]()
		- [RDP]()
 




