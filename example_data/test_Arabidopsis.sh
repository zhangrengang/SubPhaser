prefix=Arabidopsis_suecica
# download genome
url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/202/805/GCA_019202805.1_ASM1920280v1/GCA_019202805.1_ASM1920280v1_genomic.fna.gz
[ ! -s ${prefix}_genome.fasta.gz ] && \
	wget $url -O ${prefix}_genomic.fna.gz -c && \
	mv ${prefix}_genomic.fna.gz ${prefix}_genome.fasta.gz

DT=`date +"%y%m%d%H%M"`
# run subphaser
options="-pre ${prefix}_" # to avoid conficts
subphaser -i ${prefix}_genome.fasta.gz -c ${prefix}_sg.config $options 2>&1 | tee ${prefix}.log.$DT

	
