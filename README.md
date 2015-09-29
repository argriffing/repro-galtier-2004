
<!-- this table of contents is written by hand... -->
# Repro-Galtier-2004
 * [Introduction](#introduction)
 * [NHML data and code](#data-and-code)
   - [Data source](#data-source)
   - [Code source](#code-source)
 * [Data pre-processing](#data-pre-processing)
   - [Alignment extraction](#alignment-extraction)
   - [Tree estimation](#tree-estimation)
   - [Tree annotation removal](#tree-annotation-removal)
 * [NHML analysis](#nhml-analysis)
   - [Inputs](#nhml-inputs)
   - [Outputs](#nhml-outputs)
   - [Command](#nhml-command)
 * [Reimplementation analysis](#reimplementation-analysis)
   - [Inputs](#reimplementation-inputs)
   - [Outputs](#reimplementation-outputs)
   - [Command](#reimplementation-command)


<a name="introduction"/>
# Introduction

This github repo documents an attempt to reproduce parts of
the results of the following publication:

```
Markov-Modulated Markov Chains and the Covarion Process of Molecular Evolution
By N. Galtier and A. Jean-Marie.
Journal of Computational Biology. 2004, 11(4): 727-733.
Published in Volume: 11 Issue 4: January 20, 2005
doi:10.1089/cmb.2004.11.727.
http://online.liebertpub.com/doi/pdf/10.1089/cmb.2004.11.727
```



<a name="data-and-code"/>
# Links to data and code at P&ocirc;le BioInformatique Lyonnais

<a name="data-source"/>
## data source
The upstream data is a large-subunit ribosomal dna
[alignment](http://pbil.univ-lyon1.fr/datasets/gcanc/lsu_align.html).
It is an html file with asterisks indicating 'good' alignment columns.
```
$ md5sum lsu_align.html 
31be444022f40d51a1e9625d38d0c579  lsu_align.html
```

<a name="code-source"/>
## code source
The C code associated with the publication is available at
`ftp://pbil.univ-lyon1.fr/pub/mol_phylogeny/nhml`.
This is called nhml (Non-Homogeneous Maximum Likelihood?)
and is for 'maximum likelihood phylogenetic inferences from DNA
sequence data using a non-homogeneous Markov model of DNA sequence evolution.'
```
$ md5sum nhml3.tar 
8909fece00f4af6806ae219b2fc9b6a4  nhml3.tar
```


# Data pre-processing

<a name="alignment-extraction"/>
## 1) Extraction of alignment from html

todo : custom python hack was written, to make 'phylip' and 'mase' alignments

<a name="tree-estimation"/>
## 2) Rough estimation of tree from alignment

Use the [phyml online form](http://www.atgc-montpellier.fr/phyml/)
to roughly construct a phylogenetic tree for which the branch length
units are nucleotide substitution events per site.
The input file is a phylip alignment which has been pre-processed
to include only the relevant columns.
The default phyml settings on the web form are used,
except that the phylip sequence file format choice is changed from
'interleaved' to 'sequential' as appropriate.

<a name="tree-annotation-removal"/>
## 3) Removal of extraneous annotation in the phyml tree output

todo : dendropy was used


<a name="nhml-analysis"/>
# nhml analysis

<a name="nhml-inputs"/>
## inputs
 * `galtier2004.mase` (nucleotide alignment in
   [mase](http://pbil.univ-lyon1.fr/help/formats.html) format)
 * `galtier2004.tree` (phylogenetic tree in newick format
   with branch lengths but without any other annotation)
 * `galtier2004.opt` (options specific to the nhml3 software)

<a name="nhml-outputs"/>
## outputs
 * `stdout` (parameter estimates and log likelihood)
 * `treefile.eqgc` (reports branch-specific G+C process parameter estimates;
   this output is ignored)
 * `treefile.ndgc` (reports conditional G+C proportions per node;
   this output is also ignored)

<a name="nhml-command"/>
## command
```
$ /path/to/nhml3/sources/eval_nh galtier2004.mase galtier2004.tree galtier2004.opt
```


<a name="reimplementation-analysis"/>
# reimplementation analysis

<a name="reimplementation-inputs"/>
## inputs

<a name="reimplementation-outputs"/>
## outputs

<a name="reimplementation-command"/>
## command

