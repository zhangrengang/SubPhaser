prefix=peanut
# download genome
url=https://peanutbase.org/data/v2/Arachis/hypogaea/genomes/Tifrunner.gnm2.J5K5/arahy.Tifrunner.gnm2.J5K5.genome_main.fna.gz
[ ! -s ${prefix}_genome.fasta.gz ] && \
	wget $url -O ${prefix}_genomic.fna.gz -c && \
	mv ${prefix}_genomic.fna.gz ${prefix}_genome.fasta.gz

# run subphaser
options="-pre ${prefix}_" # to avoid conficts
subphaser -i ${prefix}_genome.fasta.gz -c ${prefix}_sg.config $options 2>&1 | tee ${prefix}.log
	
