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


Links to data and code at P&ocirc;le BioInformatique Lyonnais
-------------------------------------------------------------

upstream data
-------------
The upstream data is a large-subunit ribosomal dna
[alignment](http://pbil.univ-lyon1.fr/datasets/gcanc/lsu_align.html).
It is an html file with asterisks indicating 'good' alignment columns.
```
$ md5sum lsu_align.html 
31be444022f40d51a1e9625d38d0c579  lsu_align.html
```

upstream code
-------------
The C code associated with the publication is available at
`ftp://pbil.univ-lyon1.fr/pub/mol_phylogeny/nhml`.
This is called nhml (Non-Homogeneous Maximum Likelihood?)
and is for 'maximum likelihood phylogenetic inferences from DNA
sequence data using a non-homogeneous Markov model of DNA sequence evolution.'
```
$ md5sum nhml3.tar 
8909fece00f4af6806ae219b2fc9b6a4  nhml3.tar
```
