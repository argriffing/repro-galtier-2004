# http://online.liebertpub.com/doi/pdf/10.1089/cmb.2004.11.727

# NOTE: customize this path to the nhml executable.
# This code could be available on github,
# and could also be available at a lab ftp site.
# The version on github should have a more conventional termination status.
# ftp://pbil.univ-lyon1.fr/pub/mol_phylogeny/nhml
#
NHML = $$HOME/git-repos/nhml3-unofficial/nhml3/sources/eval_nh
NHML_MASE_TREE = aln.mase aln.phylip_phyml_tree.txt

# Following Table 2,
# values for GAMMA_NBCLASS (number of classes for discretized gamma rates):

NUMBERS = 3 4 5 6 8 10 15 20 50


all: mle.txt


# This first phase formats the alignment and tree.

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



# This second phase runs the analyses with the tree and alignment.

numbers.txt:
	python create-nhml-opt.py $(NUMBERS)
	touch numbers.txt

foo:
	python run-nhml.py \
		--executable=$(NHML) \
		--alignment=aln.mase \
		--tree=aln.phylip_phyml_tree.txt \
		nhml-options-*.opt


clean:
	rm --force mle.txt
