#!/usr/bin/bash
#
# This is an example of how to use DLCpar to reconcile gene families.
#
#
# Don't execute this script all at once.  Instead try copying and pasting
# each command to the command line one at a time in order to learn how the
# commands work.
#

#=============================================================================
# setup/install

# Make sure tools are compiled and installed before running the commands in
# this tutorial.  See INSTALL.txt for more information.

# Or you can run from the source directory:

cd ..
python setup.py install

cd examples
export PATH=$PATH:../bin
export PYTHONPATH=$PYTHONPATH:../python


#=============================================================================
# reconcile gene tree using DLCpar

# show help information
dlcpar -h

# Usage: dlcpar [options] <gene tree> ...
#
# Options:
#   Input/Output:
#     -s <species tree>, --stree=<species tree>
#                         species tree file in newick format
#     -S <species map>, --smap=<species map>
#                         gene to species map
#
#   File Extensions:
#     -I <input file extension>, --inputext=<input file extension>
#                         input file extension (default: "")
#     -O <output file extension>, --outputext=<output file extension>
#                         output file extension (default: ".dlcpar")
#
#   Miscellaneous:
#     -x <random seed>, --seed=<random seed>
#                         random number seed

# by default, dlcpar outputs the reconciliation in LCT format
# this creates the files 0.dlcpar{.info,.tree,.recon,.order}
dlcpar \
    -s config/paper.stree \
    -S config/paper.smap \
    -I .coal.tree -O .dlcpar \
    -x1234 \
    data/paper/0/0.coal.tree

# convert to 3T format (alternatively, use '--output_format=dlcoal' when running dlcpar)
# this creates the files 0.dlcpar{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}
# duplications and losses should be inferred using the locus tree and locus reconciliation
dlcpar_to_dlcoal -s config/paper.stree data/paper/0/0.dlcpar.tree

# clean up
find data/paper -name '*dlcpar*' | xargs rm


#=============================================================================
# DLCpar parameters

# The default costs are D=1, L=1, C=0.5.  To change these costs:
#     '-D <dup cost> -L <loss cost> -C <coal cost>'
# If you use different costs, you should also change the prescreen parameters:
#     '--prescreen_min <prescreen_min> --prescreen_factor <prescreen factor>'


# The default is to use the LCT model and bound the search space using heuristics
# (DLCpar-bound in the paper).  DLCpar as presented in the paper (with full search)
# uses the options below.
#
# (1) To run DLCpar without bounding the number of dups/losses per species branch:
#     '--max_dups=-1 --max_losses=-1'
#
# (2) To run DLCpar without prescreening the reconciliations:
#     '--no_prescreen'


#=============================================================================
# run DLCpar using search

# For some gene trees that are very large or highly incongruent to the species tree,
# you may need to use DLCpar-search, which relies on the 3T model and searches the space
# of reconciliations using a hill-climbing approach.

# show help information
dlcpar_search -h

# Usage: dlcpar_search [options] <gene tree> ...
#
# Options:
#   Input/Output:
#     -s <species tree>, --stree=<species tree>
#                         species tree file in newick format
#     -S <species map>, --smap=<species map>
#                         gene to species map
#
#   File Extensions:
#     -I <input file extension>, --inputext=<input file extension>
#                         input file extension (default: "")
#     -O <output file extension>, --outputext=<output file extension>
#                         output file extension (default: ".dlcpar")
#
#   Search:
#     -i <# iterations>, --iter=<# iterations>
#                         number of search iterations (default: 10)
#     --nprescreen=<# prescreens>
#                         number of prescreening iterations (default: 20)
#
#   Miscellaneous:
#     -x <random seed>, --seed=<random seed>
#                         random number seed

# by default, dlcpar_search outputs the reconciliation in 3T format
# this creates the files 0.dlcpar{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}
# however, this CANNOT be converted to LCT format because the locus tree is undated
# duplications and losses should be inferred using the locus tree and locus reconciliation
dlcpar_search \
    -s config/paper.stree \
    -S config/paper.smap \
    -I .coal.tree -O .dlcpar \
    -i 1000 --nprescreen 20 \
    -x1234 \
    data/paper/0/0.coal.tree

# clean up
find data/paper -name '*dlcpar*' | xargs rm


#=============================================================================
# convert between labeled coalescent tree and three-tree formats
# and infer events based on two formats

# convert from 3T to LCT
# let the input file be named <base><inputext>
# then new files are named <base><outputext> with LCT extensions
# here, this creates the files 0.dlcpar{.info,.tree,.recon,.order}
dlcoal_to_dlcpar -s config/paper.stree -S config/paper.smap data/paper/0/0.coal.tree

# convert from LCT to 3T
# let the input file be named <base><inputext>
# then new files are named <base><outputext> with 3T extensions
# here, this creates the files 0.dlcpar{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}
dlcpar_to_dlcoal -s config/paper.stree data/paper/0/0.dlcpar.tree

# find true events using 3T format
# the extension tells the script how to find the necessary files
# which need to be named <base>{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}
echo data/paper/0/0.coal.tree | tree-events-dlc -s config/paper.stree -S config/paper.smap -T .coal.tree

# find true events using LCT format
# the extension tells the script how to find the necessary files
# which need to be named <base>{.tree,.recon,.order}
echo data/paper/0/0.dlcpar.tree | tree-events-dlcpar -s config/paper.stree -S config/paper.smap -T .tree

# clean up
find data/paper -name '*dlcpar*' | xargs rm


#=============================================================================
# count number of optimal reconciliations and uniformly sample optimal reconciliations

# By default, dlcpar returns a single (uniformly sampled) optimal reconciliation.
# It also counts the number of optimal reconciliations and outputs this count
# (in <base>.dlcpar.info). To sample multiple optima, use '-n <number of reconciliations>'.
# In this case, solutions will be separated by '# Solution 0', '# Solution 1', etc.

# Note that dlcpar_search cannot sample multiple optima. Furthermore, other scripts
# (dlcpar_to_dlcoal, dlcpar_to_dlcoal, tree-events-dlc, and tree-events-dlcpar)
# expect a single solution per file, so you will have to separate the solutions
# into individual files if you use these scripts.


#=============================================================================
# find the cost space and events using DLCscape

# show help information
dlcscape -h

# Usage: dlcscape [options] <gene tree> ...
#
# Options:
#   Input/Output:
#     -s <species tree>, --stree=<species tree>
#                         species tree file in newick format
#     -S <species map>, --smap=<species map>
#                         gene to species map
#     --events {U, I}	Report (U)nion or (I)ntersection of events - if --events flag is not present,
#                       does not compute events at all, which is faster if you just want the cost space
#     --draw_regions    Draw regions to the screen? Need x-server forwarding if you are running it on a server
#                       (use ssh -Y)
#
#   File Extensions:
#     -I <input file extension>, --inputext=<input file extension>
#                         input file extension (default: "")
#     -O <output file extension>, --outputext=<output file extension>
#                         output file extension (default: ".dlcpar")
#
#   Costs:
#     -D <min>-<max>, --duprange=<min>-<max>
#                      range of duplication cost, normalized to coalescence cost = 1
#     -L <min>-<max>, --lossrange=<min>-<max>
#                      range of loss cost, normalized to coalescence cost = 1
#
#   Miscellaneous:
#     -x <random seed>, --seed=<random seed>
#                         random number seed

# this creates the files 0.dlcscape{.info,.regions,.events}
# uses the default dup/loss range of 0.2-5
dlcscape \
    -s config/paper.stree \
    -S config/paper.smap \
    -I .coal.tree -O .dlcscape \
    -x1234 --events I \
    data/paper/0/0.coal.tree

# view landscape, this creates the file 0.dlcscape.pdf
view_dlcscape -o data/paper/0/0.dlcscape.pdf data/paper/0/0.dlcscape.regions

# clean up
find data/paper -name '*dlcscape*' | xargs rm
