DLCpar
http://compbio.mit.edu/dlcpar/
Yi-Chieh Wu, with libraries contributed from Matthew Rasmussen

=============================================================================
ABOUT

DLCpar is a reconciliation program that maps a gene tree to a species tree
by inferring gene duplications, losses, and coalescence (accounting for
incomplete lineage sorting).  DLCpar uses the labeled coalescent tree (LCT)
to infer the species and locus to which a gene belongs.

DLCpar citation:
Wu, Rasmussen, Bansal, Kellis. Most Parsimonious Reconciliation
in the Presence of Gene Duplication, Loss, and Deep Coalescence
using Labeled Coalescent Trees.
Genome Research (24)3:475-486, 2014.

This package includes the Python source code of the DLCpar program,
as well as several useful utilities for working with LCTs.


=============================================================================
USAGE

Running dlcpar with no arguments will print out its command-line usage.


#=============================================================================
# File formats

Labeled Coalescent Tree (LCT) -- used by DLCpar
  X.tree        -- coalescent tree (newick)
  X.recon       -- reconciliation from coal tree to species tree and locus set
                     col1 = coal node ID
                     col2 = species node ID
                     col3 = locus
  X.order       -- partial order of coal tree nodes
                     col1 = species node ID
                     col2 = parent locus
                     col3 = (ordered) list of coal node IDs (comma-separated)


Three Tree Model (3T) -- used by DLCoalRecon
  X.coal.tree   -- coalescent tree (newick)
  X.coal.recon  -- reconciliation from coal tree to locus tree
                     col1 = coal node ID
                     col2 = locus node ID
                     col3 = event ("none")
  X.locus.tree  -- locus tree (newick)
  X.locus.recon -- reconciliation from locus tree to species tree
                     col1 = locus node ID
                     col2 = species node ID
                     col3 = event ("gene", "spec", "dup")
  X.daughters   -- daughter nodes in locus tree (one per line)


#=============================================================================
# Examples

See examples/test.sh for an example of how to use each program in this
package.
