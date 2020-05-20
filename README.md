# Overview

DLCpar is a phylogenetic program for inferring a most parsimonious reconciliation between a gene tree and species tree under a duplication-loss-coalescence model. DLCpar is a reconciliation package that maps a gene tree to a species tree by inferring gene duplications, losses, and coalescence (accounting for incomplete lineage sorting). DLCpar uses the labeled coalescent tree (LCT) to infer the species and locus to which a gene belongs.

More detail can be found at the [project website](http://www.cs.hmc.edu/~yjw/software/dlcpar/), [manual](MANUAL.md), and [examples](EXAMPLES.md).


## Authors
- Yi-Chieh Wu
- Ross Mawhorter and Nuo Liu (cost landscapes)
- Haoxing Du and Yi Sheng Ong (multiple optimal reconciliations)
- Matthew Rasmussen (compbio libraries)


## Citations

1. Wu, Rasmussen, Bansal, Kellis.\
Most Parsimonious Reconciliation in the Presence of Gene Duplication, Loss, and Deep Coalescence using Labeled Coalescent Trees.\
[Genome Research 24(3):475-486, 2014.](http://dx.doi.org/10.1101/gr.161968.113)\
This is the recommended citation for the software.

2. Du, Ong, Knittel, Mawhorter, Liu, Gross, Tojo, Libeskind-Hadas, Wu.\
Multiple Optimal Reconciliations under the Duplication-Loss-Coalescence Model.\
[IEEE/ACM Transactions on Computational Biology and Bioinformatics. In press.](http://dx.doi.org/10.1109/TCBB.2019.2922337)\
This citation should be used if you are counting the number of optimal reconciliations or using uniform random sampling.

3. Mawhorter, Liu, Libeskind-Hadas, Wu.\
Infering Pareto-Optimal Reconciliations across Multiple Event Costs under the Duplication-Loss-Coalescence Model.\
[BMC Bioinformatics 20:639, 2019.](http://dx.doi.org/10.1186/s12859-019-3206-6)\
This citation should be used if you are partitioning the landscape of event costs.


# License

Copyright (c) 2012-2019 by Yi-Chieh Wu, with modified Python libraries, original (C) 2005-2011 by Matthew Rasmussen.

DLCpar is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
