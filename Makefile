
all: aln.phylip_phyml_tree.txt

lsu_align.html:
	curl -O http://pbil.univ-lyon1.fr/datasets/gcanc/lsu_align.html

lsu_align.trimmed: lsu_align.html
	tail -n +16 lsu_align.html | head -n -3 > lsu_align.trimmed

aln.mase: lsu_align.trimmed
	python process_1.py --mase-out=aln.mase < lsu_align.trimmed

aln.phylip: lsu_align.trimmed
	python process_1.py --phylip-out=aln.phylip < lsu_align.trimmed

aln.phylip_phyml_tree.txt: aln.phylip
	phyml --sequential --input aln.phylip --bootstrap 0
