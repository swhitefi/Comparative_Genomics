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
