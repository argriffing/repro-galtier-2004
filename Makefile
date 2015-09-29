
# this will need to be customized...
NHML = $$HOME/packages/nhml3/nhml3/sources/eval_nh



all: mle.txt

# download the alignment as an html file
lsu_align.html:
	curl -O http://pbil.univ-lyon1.fr/datasets/gcanc/lsu_align.html

# trim the first 15 and last 3 lines from the html file
lsu_align.trimmed: lsu_align.html
	tail -n +16 lsu_align.html | head -n -3 > lsu_align.trimmed

# create a mase alignment with selected columns
aln.mase: lsu_align.trimmed
	python process_1.py --mase-out=aln.mase < lsu_align.trimmed

# create a phylip alignment with selected columns
aln.phylip: lsu_align.trimmed
	python process_1.py --phylip-out=aln.phylip < lsu_align.trimmed

# use phyml to build a tree without branch support annotation
aln.phylip_phyml_tree.txt: aln.phylip
	phyml --sequential --input aln.phylip --bootstrap 0

# search for maximum likelihood parameter estimates
mle.txt: aln.mase aln.phylip_phyml_tree.txt galtier2004.opt
	$(NHML) aln.mase aln.phylip_phyml_tree.txt galtier2004.opt > mle.txt

clean:
	rm --force mle.txt
