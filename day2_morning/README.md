# Day 2 Morning
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

On day 1 we worked through a pipeline to map short-read data to a pre-existing assembly and identify single-nucleotide variants (SNVs) and small insertion/deletions. However, what this sort of analysis misses is the existence of sequence that is not present in your reference. Today we will tackle this issue by assembling our short reads into larger sequences, which we will then analyze to characterize the functions unique to our sequenced genome.   

Execute the following command to copy files for this morning’s exercises to your scratch directory: 

```
cd /scratch/micro612w16_fluxod/username
cp -r /scratch/micro612w16_fluxod/shared/data/day2_morn ./
```

## Genome Assembly using Spades Pipeline
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

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
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

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
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

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
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

## Genome Annotation
[[back to top]](https://github.com/alipirani88/Comparative_Genomics#bacterial-comparative-genomics-workshop)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

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
