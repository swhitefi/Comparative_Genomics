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

## Perform QC on fastq files
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

## Examine results of SPANDx pipeline
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

## Recombination detection and tree generation
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

## Phylogenetic tree annotation and visualization
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

## Assessment of genomic deletions
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
 
=======
>>>>>>> 3cd338e0b0304d558fa4bd2935552735707d5465

```
blastall -p blastp -i Abau_all.pfasta -d ardb_ROI_genes.pfasta -o bl_blastp_results -m 8 -e 1e-20 -v 1 -b 1
```

## Identification of antibiotic resistance genes with LS-BSR and the ARDB database
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

## Perform pan-genome analysis with LS-BSR
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

** Perform genome comparisons with ACT**
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

## Perform whole genome alignment with Mauve and convert alignment to other useful formats
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

## Perform some DNA sequence comparisons and phylogenetic analysis in ape
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

## Perform SNP density analysis to discern evidence of recombination
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

## Perform recombination filtering with gubbins
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

## Perform QC on fastq files
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

## Examine results of SPANDx pipeline
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

## Recombination detection and tree generation
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

## Phylogenetic tree annotation and visualization
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

## Assessment of genomic deletions
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
