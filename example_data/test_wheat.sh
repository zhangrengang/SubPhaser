prefix=wheat
# download genome
[ ! -s ${prefix}_genome.fasta.gz ] && wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v2.1/iwgsc_refseqv2.1_assembly.fa.zip -c && \
	unzip iwgsc_refseqv2.1_assembly.fa.zip && \
	mv iwgsc_refseqv2.1_assembly.fa ${prefix}_genome.fasta && \
	gzip ${prefix}_genome.fasta -f && \
	rm iwgsc_refseqv2.1_assembly.fa.zip

# run subphaser
options="-pre ${prefix}_" # to avoid conficts
subphaser -i ${prefix}_genome.fasta.gz -c ${prefix}_sg.config $options
