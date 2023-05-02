# Processing of repeatmasker tracks

- repeatmasker tracks were downloaded from [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables) in bed format.

- the bedfiles were uncompressed and sorted:
	- `zcat repeatmasker.hg38.bed.gz | sort -k1,1 -k2,2n > repeatmasker.hg38.sorted.bed`
- the sorted bedfiles were merged:
	- `bedtools merge -i repeatmasker.hg38.sorted.bed > repeatmasker.hg38.merged.bed`
- ultimately the files were compressed via gzip
