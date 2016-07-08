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
# this creates the files 0.dlcpar{.tree,.recon,.order}
dlcpar \
    -s config/flies.stree \
    -S config/flies.smap \
    -I .coal.tree -O .dlcpar \
    -x1234 \
    sim-flies/0/0.coal.tree

# convert to 3T format (alternatively, use '--output_format=dlcoal' when running dlcpar)
# this creates the files 0.dlcpar{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}
# duplications and losses should be inferred using the locus tree and locus reconciliation
dlcpar_to_dlcoal -s config/flies.stree sim-flies/0/0.dlcpar.tree

# clean up
find sim-flies -name '*dlcpar*' | xargs rm


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
    -s config/flies.stree \
    -S config/flies.smap \
    -I .coal.tree -O .dlcpar \
    -i 1000 --nprescreen 20 \
    -x1234 \
    sim-flies/0/0.coal.tree

# clean up
find sim-flies -name '*dlcpar*' | xargs rm


#=============================================================================
# convert between labeled coalescent tree and three-tree formats

# convert from 3T to LCT
# this is only possible if the 3T has fully dated (coalescent and locus) trees
# let the input file be named <base><inputext>
# then new files are named <base><outputext> with LCT extensions
# here, this creates the files 0.dlcpar{,.tree,.recon,.order}
dlcoal_to_dlcpar -s config/flies.stree -S config/flies.smap sim-flies/0/0.coal.tree

# convert from LCT to 3T
# this is NON-REVERSIBLE because the 3T produced has non-dated (coalescent and locus) trees
# let the input file be named <base><inputext>
# then new files are named <base><outputext> with 3T extensions
# here, this creates the files 0.dlcpar{,.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}
dlcpar_to_dlcoal -s config/flies.stree sim-flies/0/0.dlcpar.tree

# clean up
find sim-flies -name '*dlcpar*' | xargs rm
