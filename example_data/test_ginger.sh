prefix=ginger
# download genome
url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/446/385/GCA_018446385.1_Zo_v1.1/GCA_018446385.1_Zo_v1.1_genomic.fna.gz
[ ! -s ${prefix}_genome.fasta.gz ] && \
	wget $url -O ${prefix}_genomic.fna.gz -c && \
	mv ${prefix}_genomic.fna.gz ${prefix}_genome.fasta.gz

# run subphaser
options="-pre ${prefix}_" # to avoid conficts
../subphaser.py -i ${prefix}_genome.fasta.gz -c ${prefix}_sg.config $options

	
