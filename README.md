
### Quick install and start ###
```
git clone https://github.com/zhangrengang/SubPhaser
cd SubPhaser

# install
conda env create -f SubPhaser.yaml
conda activate SubPhaser
python setup.py install

# start
cd example_data
# small genome (Arabidopsis_suecica: 270Mb)
bash test_Arabidopsis.sh
# middle genome	(peanut: 2.6Gb)
bash test_peanut.sh
# large genome (wheat: 14Gb)
bash test_wheat.sh
```

### Introduction ###


### Inputs ###
1. Chromosome-level genome sequences (fasta format), e.g. [the wheat genome](https://wheat-urgi.versailles.inra.fr/Seq-Repository/Assemblies) (ABD, 1n=3x=21).
2. Configuration of homologous chromosome sets, e.g. 
```
Chr1A   Chr1B   Chr1D					# each row is one homologous chromosome set
Chr2B   Chr2A   Chr2D					# seperate with blank chracter(s)
Chr3D   Chr3B   Chr3A
Chr4A   Chr4B   Chr4D
5A|Chr5A   5B|Chr5B   5D|Chr5D			# will rename chromosome id as 5A, 5B and 5D
Chr6A,Chr7A   Chr6B,Chr7B   Chr6D,Chr7D	# treat multiple chromosomes together using ","
```
### Run SubPhaser ###
Run with default parameters:
```
subphaser -i genome.fasta.gz -c sg.config
```
Run with only core algorithm enabled:
```
subphaser -i genome.fasta.gz -c sg.config -disable_ltr -disable_circos
```
### Outputs ###
```
phase-results/
├── k15_q200_f2.circos/				# config and data files for circos plot, developer may re-plot with some custom modification
├── k15_q200_f2.kmer_freq.pdf		# histogram of differential kmers, useful to adjust option `-q`
├── k15_q200_f2.kmer.mat			# differential kmer matrix (m kmer x n chromosome)
├── k15_q200_f2.kmer.mat.pdf		# heatmap of the kmer matrix
├── k15_q200_f2.kmer.mat.R			# R script for the heatmap plot
├── k15_q200_f2.kmer_pca.pdf		# PCA plot of the kmer matrix
├── k15_q200_f2.chrom-subgenome.tsv	# subgenome assignments and bootstrap values
├── k15_q200_f2.sig.kmer-subgenome.tsv	# subgenome-specific kmers
├── k15_q200_f2.bin.enrich			# subgenome-specific enrichments by genome window/bin
├── k15_q200_f2.ltr.enrich			# subgenome-specific LTR-RTs
├── k15_q200_f2.ltr.insert.pdf		# density plot of insertion age of subgenome-specific LTR-RTs
├── k15_q200_f2.ltr.insert.R		# R script for the density plot
├── k15_q200_f2.LTR_Copia.tree.pdf	# phylogenetic tree plot of subgenome-specific LTR/Copia elements
├── k15_q200_f2.LTR_Copia.tree.R	# R script for the LTR/Copia tree plot
├── k15_q200_f2.LTR_Gypsy.tree.pdf	# phylogenetic tree plot of subgenome-specific LTR/Gypsy elements
├── k15_q200_f2.LTR_Gypsy.tree.R	# R script for the LTR/Gypsy tree plot
├── k15_q200_f2.circos.pdf			# final circos plot
├── k15_q200_f2.circos.png

tmp/
├── LTR.scn                 # identifation of LTR-RTs by LTRhavest and/or LTRfinder
├── LTR.inner.fa            # inner sequences of LTR-RTs
├── LTR.inner.fa.cls.*      # classfication of LTR-RTs by TEsorter
├── LTR.filtered.LTR.fa     # full sequences of filtered LTR-RTs
├── LTR.LTR_*.aln           # alignments of LTR-RTs' protein domains
├── LTR.LTR_*.rooted.tre    # phylogenetic tree files
├── LTR.LTR_*.map           # information for tips on the tree
.....
```
An example of output figures of wheat (ABD, 1n=3x=21):
![wheat](example_data/wheat_figures.png)
