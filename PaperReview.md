# STAT 540: Paper Review Assignment
### Matthew Shannon
#### February 20th, 2019
---
For this assignment, the paper: [Gene Expression by Mouse Inner Ear Hair Cells during Development](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4405555/pdf/zns6366.pdf) published in the *Journal of Neuroscience* by Scheffer et al. in 2015 will be reviewed.

## The Goal of this Study:
This study sought to understand the gene expression of mouse Hair Cells (HCs) during development. HCs are terminally differentiated mechanoreceptor cells necessary for hearing and balance. HCs are difficult to study and defining how development influences their gene expression can greatly enhance our understanding of HC development, hereditary deafness (HD), and possible avenues for HD genetic therapies. Here, Scheffer et al. used FACS to isolate HCs from a mouse population expressing enhanced green fluorescent protein driven by the HC specific *Pou4f3* promoter. HCs were collected from the cochlear and utricle (vestibluar) tissue of male and female mice at four developmental time points occurring before and during the acquisition of mechanosensitivity: E16, P0, P4, and P7 and compared against supporting cells (SCs) of the sensory epithelium and nearby inner ear epithelial cells. Imemdiately following isolation, HCs and SCs were subjected to next-generation high-throughput sequencing (HTS) analysis to quantify transcriptomic differences present in the cells of each condition. Differnetial gene expression analysis was conducted on the RNA-seq data to compare gene expression in HC and SC cells across developmental stages followed by a principle component analysis, and a gene ontology analysis of the biological processes assoicated with HC development. From this, Scheffer et al. identified that; at E16, HCs express genes related to morphogenesis, differentiation, and development; at P0 and P4, HCs express specifc genes related to sensory and inner ear differentiation, and; at P7 HCs express genes related to cell signalling and synapsis, all relative to the gene expression profiles of SCs. 

RNA-seq findings were reproduced in a second mouse strain expressing the fluorescent tdTomatoe protein under the control of the *Gfi1* promoter, and validated using quantitative real-time polymerase chain reaction to confirm differential gene expression, as well as using *in situ* immunohistochemistry and immunocytochemistry analyses. 

From this data, Scheffer et al. concluded that genes preferentially expressed in HCs are good candidates for unknown deafness genes. This conclusion is supported by the previous identification of new deafness genes by researchers using their published data in the SHIELD database, such as the identification of the Usher syndrome gene *CIB2* by Riazuddin et al., in 2012.

## Datasets used in this Study:
Raw data from this study are available on the National Centre for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) database repository under the accession number: [GSE60019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60019). In addition,to make their data publically available, Scheffer et al. created the Shared Inner Ear Laboratory Database ([SHIELD](https://shield.hms.harvard.edu)) database which presents gene expression data integrated with compregensive annotation.

In this study, mRNA was collected from HC and SC cells of the cochlear and utricle (vestibluar) tissue of male and female mice at four developmental time points: E16, P0, P4, and P7. mRNA was then amplified and converted into a nondirectional Illumina sequencing library for single-end sequencing of 3' tagged 35/50bp reads via Illumina HiSeq methods. A comprehensive, cell-type specific RNA-seq analysis was done on the 16 samples using a matrix designed of three factors: [1] *cell type*, [2] *tissue source*, and [3] *developmental stage*. In this study biological duplicates were not used, rather, samples of different levels in other factors were combined and treated as replicates for statistical analysis.

## Analyses used in this Study:


## Paper Review:
