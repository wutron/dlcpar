original Algorithms:

DLCpar
http://www.cs.hmc.edu/~yjw/software/dlcpar/

xscape
https://www.cs.hmc.edu/xscape/

Authors
Ran Libeskind-Hadas, Jessica Yi-Chieh Wu, Mukul Bansal (xscape)
Yi-Chieh Wu(DLCpar)
Matthew Rasmussen (compbio libraries)
Ross Mawhorter and Ivy Liu (DLCscape)

=============================================================================
ABOUT

DLCscape is a phylogenetic tree reconciliation program that maps a gene tree to a species tree
by inferring gene duplications, losses, and coalescence (accounting for
incomplete lineage sorting).Given a gene tree, species tree and a mapping
from gene to species, DLCscape uses the labeled coalescent tree (LCT)
to infer the species and locus to which a gene belongs. DLCscape finds optimal
solutions without specific costs for the three evolutionary events.
With given range of costs for events, it divides the space of costs
into regions, each corredponding to some optmal reconciliations.

DLCscape present the individual events in the reconciliations in each
regions and how many reconciliations they appear in. DLCscape analyzes the
number of regions each event appears in to give event support. It also provide
users with information about the regions cost space divides into.

This package includes the Python source code of the DLCscape program
(dlcscape.py),
as well as several useful utilities for working with LCTs.

=============================================================================
CITATION

=============================================================================
LIBRARIES

rasmus 
compbio
dlcpar

=============================================================================
USAGE

Running 'dlcpar -h' will print out its command-line usage.

Options:
   Input/Output:
    -s <species tree>, --stree=<species tree> 
                         species tree file in newick format
    -S <species map>,  --smap=<species map>
                         gene to species map
    --lmap=<locus map> gene to locus map(species specific)
    --events=<output events> 
#=============================================================================
# File formats

X.tree        -- gene tree file (newick)
x.stree       -- species tree file (newick)
x.smap        -- species map
                   col1=gene ID
		   col2=species ID

#=============================================================================
# Examples

See examples/test.sh for an example of how to use each program in this
package.
