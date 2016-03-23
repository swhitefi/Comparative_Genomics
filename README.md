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
	- [Determine which genomes contain beta-lactamase genes](https://github.com/alipirani88/Comparative_Genomics#getting-your-data-onto-glux-and-setting-up-environment-variable)
	- [Identification of antibiotic resistance genes with LS-BSR and the ARDB database](https://github.com/alipirani88/Comparative_Genomics#quality-control-using-fastqc)
	- [Perform pan-genome analysis with LS-BSR](https://github.com/alipirani88/Comparative_Genomics#quality-trimming-using-trimmomatic)
	- [Perform genome comparisons with ACT](https://github.com/alipirani88/Comparative_Genomics)










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
> FastQC can be run in two modes: "command line" or as a GUI (graphical user interface). We will be using command line version of it.

>iv. Get an interactive cluster node to start running programs

>v. Run FastQC to generate quality report of sequence reads.

```
fastqc -o Rush_KPC_266_FastQC_results/before_trimmomatic/ Rush_KPC_266_1_combine.fastq.gz Rush_KPC_266_2_combine.fastq.gz --extract
```

This will generate the results directory for forward and reverse fastq reads called Rush_KPC_266_1_combine_fastqc and Rush_KPC_266_2_combine_fastqc in output folder provided with -o argument. The summary.txt file in these directories indicates if the data passed different quality control tests. You can visualize and assess the quality of data by opening html report in a local browser.


>vi. Exit your cluster node so you don’t waste cluster resources and $$$!

>vii. Download FastQC report to your home computer to examine

```
scp username@flux-xfer.engin.umich.edu:/scratch/micro612w16_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/before_trimmomatic/*.html /path-to-local-directory/
```

The analysis in FastQC is performed by a series of analysis modules. The left hand side of the main interactive display or the top of the HTML report show a summary of the modules which were run, and a quick evaluation of whether the results of the module seem entirely normal (green tick), slightly abnormal (orange triangle) or very unusual (red cross). 

`Screenshots explanation.`
`Explaining Summary results, Basic statistics, per base sequence quality, overrepresented sequences(adapters) from before trimmomatic report.`

Notice the quality drop(per base sequence quality graph) at the end of Rush_KPC_266_2_combine_fastqc.html report. This is commonly observed in illumina samples that as the number of sequencing cycles performed is increased the average quality of the base calls, as reported by the Phred Scores produced by the sequencer falls. Check the overrepresented sequences graph and the kind of adapters that were used for sequencing these samples.

Check out [this](https://sequencing.qcfail.com/articles/loss-of-base-call-accuracy-with-increasing-sequencing-cycles/) for more detailed explaination as to why quality drops with increasing sequencing cycles.

> [A video FastQC walkthrough created by FastQC developers](https://www.youtube.com/watch?v=bz93ReOv87Y "FastQC video") 

## Quality Trimming using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic "Trimmomatic Homepage")
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)

Now we will run Trimmomatic on these raw data to remove low quality reads and adapters and run FastQC again to check if the quality improved after running trimmomatic. 

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

explaining parameters and its default value. Adapter file. Changing only SLIDINGWINDOW parameter from default 4:15 to 4:20 for raw reads.

>iv. Run the below trimmomatic commands on raw reads(explaining parameters and its default value. Adapter file.)

```
time java -jar /scratch/micro612w16_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_266_1_combine.fastq.gz Rush_KPC_266_2_combine.fastq.gz Rush_KPC_266_trimmomatic_results/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results/forward_unpaired.fq.gz Rush_KPC_266_trimmomatic_results/reverse_paired.fq.gz Rush_KPC_266_trimmomatic_results/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w16_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:40 HEADCROP:0
```
We are using the default parameters changing only SLIDINGWINDOW from 4:15 to 4:20.

>v. Now create new directories in day1_morn folder and Run FastQC on these trimmomatic results.

```
mkdir Rush_KPC_266_FastQC_results/after_trimmomatic

fastqc -o Rush_KPC_266_FastQC_results/after_trimmomatic/ Rush_KPC_266_trimmomatic_results/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results/reverse_paired.fq.gz --extract
scp username@flux-xfer.engin.umich.edu:/scratch/micro612w16_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/after_trimmomatic/*.html /path-to-local-directory/
```

`screenshots explaination`
`before trimmomatic and after trimmomatic explanation: How quality and overrepresented sequences red cross signal disappeared after running trimmomatic`

The head bases in per base sequence content graph are slightly imbalanced. This is not very bad but you can fix this by running trimming these imbalanced head bases using HEADCROP:9 parameter in the above command.

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

**1. Copy day1_after directory from shared data directory in your home directory.

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

Lets try to get some statistics about various outputs that were created using these pipeline steps and check if everything makes sense.

>i. Reads Alignment statistics:
 
```
samtools flagstat Rush_KPC_266__aln.bam > Rush_KPC_266__alignment_stats
```

Explain the output!

ii. VCF statistics:  

```
bgzip Rush_KPC_266__aln_mpileup_raw.vcf   
tabix Rush_KPC_266__aln_mpileup_raw.vcf.gz  
vcf-stats Rush_KPC_266__aln_mpileup_raw.vcf.gz > Rush_KPC_266__raw_vcf_stats
```

Explain the output!

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
Explain the output!

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

# Day 2 – Afternoon

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

**1. Determine which genomes contain beta-lactamase genes**

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

**2. Identification of antibiotic resistance genes with LS-BSR and the ARDB database**

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

**3. Perform pan-genome analysis with LS-BSR**

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



