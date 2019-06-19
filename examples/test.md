These are examples of how to use DLCpar to reconcile gene families.

Do not execute this script all at once. Instead, try copying and pasting each command to the command line one at a time in order to learn how thte commands work.

# Setup/install

Make sure tools are compiled and installed before running the commands in this tutorial. See INSTALL.txt for more information.

Or you can run from the source directory:

    cd ..
    python setup.py install

    cd examples
    export PATH=$PATH:../bin
    export PYTHONPATH=$PYTHONPATH:../python

# Usage

## Reconcile gene tree using DLCpar dynamic programming

### DLCpar dynamic programming

 By default, dlcpar outputs the reconciliation in LCT format. This creates the files `0.dlcpar{.info,.lct.tree,.lct.recon,.lct.order}`

    dlcpar dp \
        -s config/paper.stree \
        -S config/paper.smap \
        -I .coal.tree -O .dlcpar \
        -x1234 \
        data/paper/0/0.coal.tree

### Visualise tree
This creates a svg file based on the reconciliation from `dlcpar{.lct.tree,.lct.recon,.lct.order}`

    dlcpar view_recon -s config/paper.stree --xscale 10 data/paper/0/0.dlcpar -o data/paper/0/0.dlcpar.svg

### Convert to 3T format: 
(alternatively, use '--output_format=dlcoal' when running dlcpar)

This creates the files `0.dlcpar. {.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}`

Duplications and losses should be inferred using the locus tree and locus reconciliation.

    dlcpar convert --lct_to_3tree -s config/paper.stree -S config/paper.smap data/paper/0/0.dlcpar.lct.tree

### Clean up
    find data/paper -name '*dlcpar*' | xargs rm

The default is to use the LCT model and bound the search space using heuristics (DLCpar-bound in the paper).  

## Run DLCpar using search

For some gene trees that are very large or highly incongruent to the species tree, you may need to use DLCpar-search, which relies on the 3T model and searches the space of reconciliations using a hill-climbing approach.

### DLCpar search
By default, dlcpar_search outputs the reconciliation in 3T format.

This creates the files `0.dlcpar{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}`

Duplications and losses should be inferred using the locus tree and locus reconciliation.

    dlcpar search \
        -s config/paper.stree \
        -S config/paper.smap \
        -I .coal.tree -O .dlcpar \
        -i 1000 --nprescreen 20 \
        -x1234 \
        data/paper/0/0.coal.tree

### clean up
    find data/paper -name '*dlcpar*' | xargs rm

## Convert between labeled coalescent tree and three-tree formats

### Convert from 3T to LCT

By default, `inputext=".coal.tree"` and `outputext="".`

Here, this creates the files `0.{lct.tree,lct.recon,lct.order}`:

    dlcpar convert --3tree_to_lct -s config/paper.stree -S config/paper.smap data/paper/0/0.dlcpar.coal.tree

### Convert from LCT to 3T

By default, `inputext=".lct.tree"` and `outputext=""`.

Here, this creates the files `0.{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}`:

    dlcpar convert --lct_to_3tree -s config/paper.stree -S config/paper.smap data/paper/0/0.dlcpar.lct.tree

## Infer events based on two formats

### Find true events using LCT format

The extension tells the script how to find the necessary files which need to be named `"base"{.lct.tree,.lct.recon,.lct.order}`

    dlcpar events --format lct -s config/paper.stree -S config/paper.smap -I .lct.tree data/paper/0/0.dlcpar.lct.tree

### Find true events using 3T format

The extension tells the script how to find the necessary files which need to be named `"base".{coal.recon,locus.tree,locus.recon,daughters}`

    dlcpar events --format 3t -s config/paper.stree -S config/paper.smap -I .coal.tree data/paper/0/0.dlcpar.coal.tree

## Finding equalities of trees:

### Find equality in 3t format

Suppose we have two trees, 0.coal.tree and 0.dlcpar.coal.tree, which we want to check the equality for. If they are from the same species tree and species map, we may use the following command to check their equality. If it is true, it will return 'True'. Otherwise, it will result in a key error.

    dlcpar equal --format 3t -s config/paper.stree -S config/paper.smap data/paper/0/0 data/paper/0/0.dlcpar

### Find equality in lct format

Suppose we have two trees, 0.lct.tree and 0.dlcpar.coalctl.tree, which we want to check the equality for. If they are from the same species tree and species map, we may use the following command to check their equality. If it is true, it will return 'True'. Otherwise, it will result in a key error.

    dlcpar equal --format lct -s config/paper.stree -S config/paper.smap data/paper/0/0 data/paper/0/0.dlcpar

## Find the cost space and events using DLCscape

### Making the DLCscape cost space and events files:

    dlcpar landscape \
        -s config/paper.stree \
        -S config/paper.smap \
        -I .coal.tree -O .dlcscape \
        -x1234 --events I\
        data/paper/0/0.coal.tree    

## View landscape

This creates the file `0.dlcscape.pdf`, a visual representation of the cost space.

    view_dlcscape -o data/paper/0/0.dlcscape.pdf data/paper/0/0.dlcscape.regions