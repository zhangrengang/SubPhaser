prefix=poplar
# download genome
[ ! -s ${prefix}_genome.fasta.gz ] && \
	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/804/465/GCA_018804465.1_PTv2/GCA_018804465.1_PTv2_genomic.fna.gz -c && \
	mv GCA_018804465.1_PTv2_genomic.fna.gz ${prefix}_genome.fasta.gz

# run subphaser
options="-pre ${prefix}_" # to avoid conficts
subphaser -i ${prefix}_genome.fasta.gz -c ${prefix}_sg.config $options

	
