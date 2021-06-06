# GnRH Receptor Tree
R code and associated files for constructing a cladogram depicting the protein sequence similarity of Gonadotropin-releasing hormone receptor. The following files are required for constructing the tree:

* *akh.tree.final.R* - R code for constructing cladogram of Gonadotropin-releasing hormone receptors. 
* *20210312.ncbi-protein-'Gonadotropin-releasing hormone receptor'.fasta* - fasta file containing all ref sequences returned from [NCBI protein](https://www.ncbi.nlm.nih.gov/protein/) returned when "Gonadotropin-releasing hormone receptor" was searched for on 2021-03-12. 
* *NP_648571.1-dmCrzR.fasta* - fasta file containing ref sequence for Drosophila melanogaster's corazonin receptor. 
* *GnRHreceptors.droplist.csv* - Comma seperated value file consistating of a column of accession numbers and booleans indicating if fasta entry should be inclulded. Sequnces annotated as pseudogene, partial, or low quality were not included in the tree. Additionally, to reduce the size of the the cladogram and enchance readability, all but one isoform of receptors were removed to colpase clades that consistant entirely of receptor isoforms of from the same species. 

![Protein Sequence Similarity of Gonadotropin-Releasing Hormone Receptors](https://github.com/JakeSaunders/GnRHRtree/blob/main/GnRHReceptor.tree.jpg)
