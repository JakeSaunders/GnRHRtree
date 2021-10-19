# GnRH Receptor Tree

Featured in: [JM Nelson, CJ Saunders, EC Johnson. 2021. The intrinsic nutrient sensing adipokinetic hormone producing cells function in modulation of metabolism, activity, and stress. Int. J. Mol. Sci, 22, 7515. doi.org/10.3390/ijms22147515](doi.org/10.3390/ijms22147515)

R code and associated files for constructing a cladogram depicting the protein sequence similarity of Gonadotropin-releasing hormone receptor. The following files are required for constructing the tree:

* **akh.tree.final.R** - R code for constructing cladogram of Gonadotropin-releasing hormone receptors. 
* **20210312.ncbi-protein-'Gonadotropin-releasing hormone receptor'.fasta** - fasta file containing all ref sequences returned from [NCBI protein](https://www.ncbi.nlm.nih.gov/protein/) returned when "Gonadotropin-releasing hormone receptor" was searched for on 2021-03-12. 
* **NP_648571.1-dmCrzR.fasta** - fasta file containing ref sequence for *Drosophila melanogaster*'s corazonin receptor. 
* **GnRHreceptors.droplist.csv** - Comma seperated value file consistating of a column of accession numbers and booleans indicating if fasta entry should be inclulded. Sequnces annotated as pseudogene, partial, or low quality were not included in the tree. Additionally, to reduce the size of the the cladogram and enchance readability, all but one isoform of receptors were removed to colpase clades that consistant entirely of receptor isoforms of from the same species. 
* **key_entrez.txt** - A text file containing an api key allowing Rentrez to access NCBI, must be provide by user. See [Entrez Programming Utilities Help](https://www.ncbi.nlm.nih.gov/books/NBK25500/). 
* **species.info.RData** - RData file containing saved Rentrez output used for the published verison of this tree. 

![Protein Sequence Similarity of Gonadotropin-Releasing Hormone Receptors](https://github.com/JakeSaunders/GnRHRtree/blob/main/GnRHReceptor.tree.jpg)


# Works Cited

* G Yu. Using ggtree to visualize data on tree-like structures. Current Protocols in Bioinformatics, 2020, 69:e96. doi: 10.1002/cpbi.96.
* ES Wright. Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R. The R Journal, 2016, 8(1), 352-359.

Academic Free License ("AFL") v. 3.0 https://opensource.org/licenses/AFL-3.0
