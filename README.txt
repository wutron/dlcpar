DLCpar
http://www.cs.hmc.edu/~yjw/software/dlcpar/

Authors
Yi-Chieh Wu
Ross Mawhorter and Ivy Liu (cost landscapes)
Haoxing Du and Yi Sheng Ong (multiple optimal reconciliations)
Joseph Gardi and Fiona Plunkett (ILP formulation)
Matthew Rasmussen (compbio libraries)

=============================================================================
ABOUT

DLCpar is a reconciliation package that maps a gene tree to a species tree
by inferring gene duplications, losses, and coalescence (accounting for
incomplete lineage sorting). DLCpar uses the labeled coalescent tree (LCT)
to infer the species and locus to which a gene belongs.

The two main programs are DLCpar and DLCscape. The DLCpar program generates
the most parsimonious reconciliation given the costs for the three event types,
while the DLCscape program divides the cost space into regions. In each region,
the total event counts for each reconciliation are the same.

This package includes the Python source code of the DLCpar and DLCscape programs
as well as several useful utilities for working with LCTs.

=============================================================================
CITATION

In general, cite the following paper:
Wu, Rasmussen, Bansal, Kellis.
Most Parsimonious Reconciliation in the Presence of Gene Duplication, Loss,
and Deep Coalescence using Labeled Coalescent Trees.
Genome Research (24)3:475-486, 2014.

If you are counting the number of optimal reconciliations or using uniform
random sampling:
Du, Ong, Knittel, Mawhorter, Liu, Gross, Tojo, Libeskind-Hadas, Wu.
Multiple Optimal Reconciliations under the Duplication-Loss-Coalescence Model.
IEEE/ACM Transactions on Computational Biology and Bioinformatics. In press.

If you are computing landscapes:
Mawhorter, Liu, Libeskind-Hadas, Wu.
Infering Pareto-Optimal Reconciliations across Multiple Event Costs
under the Duplication-Loss-Coalescence Model
In prep.


If you are using the ILP formulation:
TBD

=============================================================================
USAGE

Running 'dlcpar -h' or 'dlcscape -h' will print out its command-line usage.


=============================================================================
FILE FORMATS

# File formats for reconciliations

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


# File formats for landscapes and events over landscapes
 X.regions      -- region output file
                      Each line has a shapely region object with the position
		      of the region and the area of the region. Use
		      view_regions to display the cost landscape for a region
		      file.

 X.events       -- event output file
                       First part of the file:
		        Each line starts with an event count vector,
			followed by the count of optimal reconciliations
			with that count vector and the events in all of those
			reconciliations
		       Second part of the file:
		        Each line starts with a number, followed by a list of
			events appearing in that number of regions.


# File format of logs
 X.info         -- information output file, with dlcpar or dlcscape version,
                   program runtime and command.


=============================================================================
EXAMPLES

See examples/test.sh for an example of how to use each program in package.

