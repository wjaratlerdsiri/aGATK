# Assembly-based Genome Analysis Toolkit (aGATK)

![](/data/aGATK-workflow-v2.png)

# Main steps

We begin by generating an individual’s genome draft and extracting sample-specific target sequences. Next, we generate cohort-wise multiple sequence alignment (MSA) per gene and perform subsequent local and amino acid-guided alignment for accurate estimation of sequence qualities prior to use in tests for positive selection. Assembly-based genome analysis workflow consists of the following steps:

## Perform de novo genome assembly:

For each genome, filtered FASTQ reads generated in Supplementary Information 2.3 were assembled with pair information using ABySS v2.0.2 (Simpson et al., 2009; Supplementary Information). Default paired-end sequencing parameters for human genome assemblies were set with k-mer equal to 96 and five pairs minimally required to construct a contig. Subsequent assembled scaffolds in a FASTA format were chosen if larger than one kb.

## Locally align to sample-specific references:

The sequence assembly per sample is treated as a reference genome. For each individual, nucleotide sequences of targeted regions of interest, consensus coding sequences (CCDS; release 20) described in SI10 are locally aligned in a BAM format against the sample-specific reference using the BWA-MEM algorithm, fast and accurate local alignment (Li & Durbin, 2009; Supplementary Information). The BAM alignment is then sorted and indexed.

```
bwa mem -t 16 $REFERENCE CCDSv20.GRCh38.exon.nt.ALL.fasta \
       -R '@RG\tID:CCDS.GRCh38.exon.nt\tPL:None\tPU:None\tLB:1\tSM:CCDSv20\tCN:hcpcg' \
       | samtools view -h -S -b - > $REFERENCE.CCDSv20.GRCh38.exon.bam     

samtools sort -o $REFERENCE.CCDSv20.GRCh38.exon.sort.bam $REFERENCE.CCDSv20.GRCh38.exon.bam      

samtools index $REFERENCE.CCDSv20.GRCh38.exon.sort.bam
```

## Extract sample-specific target sequences:

At the individual level, we extract genomic coordinates and their subsequent nucleotide sequences from the BAM alignment between genome assembly and exon database using bedtools v2.26.0 (Quinlan & Hall, 2010; Supplementary Information). Targeted exon sequences aligned with each genome assembly are removed if soft-clipped.

```
bedtools bamtobed -i $BAM -cigar > $BED

bedtools getfasta -fi $REFERENCE -bed $BED -name -tab -s -fo $OUT.tab
```

## Concatenate target sequences into a gene level:

For each individual, target sequences (exons) are ordered based on exon numbers and their 5' direction and concatenated into a gene level for 18,680 genes targeted per sample, depending on an individual’s genome assembly. Multiple FASTA sequences, each of which corresponds to a single gene, are generated per sample.

## Cohort-wise multiple sequence alignment per gene:

This step consists of consolidating each gene sequence across multiple samples in order to generate per-gene MSA within a cohort. Using MAFFT v7.505 and following an installation manual (Katoh, Kuma, Toh, & Miyata, 2005; Supplementary Information), this MSA per gene is locally aligned with a corresponding reference sequence from the CCDS database for quality control and protein translation purposes in the cohort analysis in the following steps. The example data can be found in the data directory (data/OR1J1_nostop.fasta) and run by replacing the MSA variable with the path to the example file. The script `aGATK_mafft.py` written in Python programming languages can be run per gene in this step after installing the mafft binary in the path. All the Python scripts can be cloned and tested with the example data, which will take a few minutes to complete. 

```
python aGATK_mafft.py --fasta $MSA --output $OUT
```

## Multiple sequence alignment and orientation quality controls:

The next step is to quality check per-gene MSA datasets based on overall nucleotide diversity (mean p-distance) and sequence orientation relative to their hg38 reference sequences, using MEGA-CC v11.0.13 (Kumar et al., 2012; Supplementary Information) and LASTZ v1.04 programs (Harris, 2007; Supplementary Information), respectively. The MEGA program may be buggy with Apple M1 chip and needs an input file in data/ directory. There are the two following Python scripts used for these analyses:

```
python aGATK_megacc.py --target $MSA --output $OUT

python aGATK_lastz.py --target $MSA --query $MSA --output $OUT --identity 90
```

For each gene of interest, this step could be followed by strand correction, removal of gaps and sequence alignment of genes less than 100 bp. We chose only genes with >90% sequence identity among individuals studied for further analyses after repeating the above-mentioned local alignment step.

## Perform amino acid-guided alignment:

To facilitate amino acid based selection tests explained in SI10, amino acid-guided alignment of each per-gene MSA is conducted using the transAlign program (Bininda-Emonds, 2005; Supplementary Information). The step would provide accurate MSA of intact open reading frames within a cohort relative to reference coding sequences, and sequence alignment based on amino acids is often superior to that obtained directly from nucleotides. Overall nucleotide diversity is recalculated, with the MSA greater than 90% identity retained.

```
perl transAlign.pl –d$MSA -if -v -pclustalw2 -ra
```

## Preparing alignment data for selection analysis:

Formatting data for selection tests generally involves slight differences from the steps described above. This consists of (i) removing reference coding sequences (CCDS release 20) and (ii) replacing ambiguous bases (IUPAC codes) and stop codons with gaps. 

## Test for positive selection:

This step is mentioned in detail in the manuscript.


# Caveats

We perform de novo genome assembly as an alternative to short read mapping implemented in the classical GATK pipeline or as an reference genome-free analysis of short read sequencing. Despite WGS coverage of ~41X assembled in this study, target sequences of interest in some individuals are challenging to obtain, for example, 745 genes partly observed among 37 genomes analysed in Supplementary Information 10.

For the current workflow of aGATK, quality scrutiny of MSA of extracted targets is done at the gene level (2,770 genes filtered out).


