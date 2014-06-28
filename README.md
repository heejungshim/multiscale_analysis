

This repository contains scripts and (limited) results from multi-scale analyses of simulated data or real data.

## Overveiw

We have multiple ongoing projects that perform analyses of simulated data or different types of high-throughput sequencing data (e.g., DNase-seq ([Boyle et al., 2008](http://www.ncbi.nlm.nih.gov/pubmed/18243105); [Hesselberth et al., 2009](http://www.ncbi.nlm.nih.gov/pubmed/19305407)), ATAC-seq ([Buenrostro et al., 2013](http://www.ncbi.nlm.nih.gov/pubmed/24097267)), Ribo-seq ([Ingolia et al., 2011](http://www.ncbi.nlm.nih.gov/pubmed/22056041)), RNA-seq ([Mortazavi et al., 2008](http://www.ncbi.nlm.nih.gov/pubmed/18516045); [Wang et al., 2008](http://www.ncbi.nlm.nih.gov/pubmed/18978772); [Marioni et al., 2008](http://www.ncbi.nlm.nih.gov/pubmed/18550803)), ChIP-seq ([Johnson et al., 2007](http://www.ncbi.nlm.nih.gov/pubmed/17540862); [Barski et al., 2007](http://www.ncbi.nlm.nih.gov/pubmed/17512414); [Mikkelsen et al., 2007](http://www.ncbi.nlm.nih.gov/pubmed/17603471)) data) by using different multi-scale approaches (e.g., [WaveQTL](https://github.com/heejungshim/WaveQTL), [multiseq](https://github.com/stephenslab/multiseq), and [WaveHMT](https://github.com/stephenslab/hmt)). Those multi-scale analyses share scripts and results to some extent, so we have tried to put them together in one repository. For now, this repository is mostly for our collaborators to share scripts/results and replicate analyses. As some data sets are not publicly available and all analyses are work in progress, we put only limited results in the repository. However, once data sets become publicly available and projects are close to be complete, we'll share all scripts/data/results. If you are interested in our analyses (results, contributing to analyses, performing similar analyses for other applications), contact hjshim@gmail.com. 

### Compare multiseq to WaveQTL on simulated data

We simulated null and alternative data sets from the 578 dsQTLs, identified by either or both of the wavelet-based or 100bp window approach at FDR=0.01 in [Shim and Stephens 2014](https://github.com/heejungshim/WaveQTL/tree/master/doc/paper), with a procedure similar to those described in [Supplementary Material of Shim and Stephens 2014](https://github.com/heejungshim/WaveQTL/tree/master/doc/paper). We then compare performances of two approaches ([WaveQTL](https://github.com/heejungshim/WaveQTL) and [multiseq](https://github.com/stephenslab/multiseq)) on simulated data sets with different sample sizes or different read depths.


