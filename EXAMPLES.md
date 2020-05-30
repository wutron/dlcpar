<h1>Table of Contents</h1>

These are examples of how to use DLCpar commands. Try copying and pasting each command to the command line one at a time to learn how the commands work.

<!-- TOC -->

- [1. Setup / Install](#1-setup--install)
- [2. Examples](#2-examples)
    - [2.1 Reconcile gene tree using dynamic programming](#21-reconcile-gene-tree-using-dynamic-programming)
        - [Infer an MPR](#infer-an-mpr)
        - [Visualize a reconciliation on screen](#visualize-a-reconciliation-on-screen)
        - [Setting commonly-used parameters](#setting-commonly-used-parameters)
    - [2.2 Reconcile gene tree using heuristic search](#22-reconcile-gene-tree-using-heuristic-search)
        - [Infer an MPR](#infer-an-mpr-1)
    - [2.3 Convert between reconciliation formats](#23-convert-between-reconciliation-formats)
        - [Convert from three-tree to LCT](#convert-from-three-tree-to-lct)
        - [Convert from LCT to three-tree](#convert-from-lct-to-three-tree)
    - [2.4 Infer events from reconciliations](#24-infer-events-from-reconciliations)
        - [Infer from LCT](#infer-from-lct)
        - [Infer from three-tree](#infer-from-three-tree)
    - [2.5 Compare reconciliations](#25-compare-reconciliations)
        - [Compare three-trees](#compare-three-trees)
        - [Compare LCTs](#compare-lcts)
    - [2.6 Find MPR landscapes](#26-find-mpr-landscapes)
        - [Find equivalent regions and events per region](#find-equivalent-regions-and-events-per-region)
        - [View landscape](#view-landscape)
    - [Clean up](#clean-up)

<!-- /TOC -->




<!-- Install -->

# 1. Setup / Install

Make sure tools are compiled and installed before running the commands in this tutorial. See the [manual](MANUAL.md#1-install) for more information.

Or you can run from the source directory by setting these environment variables:

    cd examples
    export PATH=$PATH:../bin
    export PYTHONPATH=$PYTHONPATH:../python

<!-- /Install -->





<!-- Usage -->

# 2. Examples

## 2.1 Reconcile gene tree using dynamic programming

### Infer an MPR

    dlcpar dp \
        -s config/paper.stree \
        -S config/paper.smap \
        -I .coal.tree -O .dlcdp \
        -x1234 \
        data/paper/0/0.coal.tree

This creates the files `0.dlcdp{.info,.lct.tree,.lct.recon,.lct.order}`.

Often, you might be interested in the reconciliation in three-tree format. You have two options: (1) set `--output_format 3t`, or (2) infer a reconciliation in LCT format (as above), then convert to three-tree format.

    dlcpar convert \
        --lct_to_3t \
        -s config/paper.stree \
        -S config/paper.smap \
        data/paper/0/0.dlcdp.lct.tree

This creates the files `0.dlcdp.{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}`.

### Visualize a reconciliation on screen

Currently, you can only view reconciliations in LCT format.

    dlcpar view_lct \
        -s config/paper.stree \
        --xscale 10 \
        data/paper/0/0.dlcdp.lct.tree

### Setting commonly-used parameters 

By default, `dlcpar dp` infers MPRs using default event costs of `D=1`, `L=1`, `C=0.5` and `K=C`. To change these costs (to any positive real number):

    -D <dup cost> -L <loss cost> -C <coal cost> -K <coal dup cost>

If you use different costs, you should also change the prescreen parameters:

    --prescreen_min <prescreen_min> --prescreen_factor <prescreen_factor>

For efficiency, the search space of reconciliations is restricted via heuristics (DLCpar-bound in the original paper). To remove the heuristics (and explore the full space):

    --max_dups -1 --max_losses -1 --no_prescreen

By default, `dlcpar dp` returns a single MPR (sampled uniformly at random). To sample multiple optima:

    -n <number of reconciliations>

Solutions are separated by "# Solution 0", "# Solution `", etc. Note that all utility and visualization commands expect a single solution per file, so the multiple solutions will have to be separated to use these commands.



## 2.2 Reconcile gene tree using heuristic search

For some gene trees that are very large or highly incongruent to the species tree, an exhaustive search may not be possible. Instead, you may need to use a heuristic hill-climbing search.

### Infer an MPR

    dlcpar search \
        -s config/paper.stree \
        -S config/paper.smap \
        -I .coal.tree -O .dlcsearch \
        -i 100 --nprescreen 20 \
        -x1234 \
        data/paper/0/0.coal.tree

This creates the files `0.dlcsearch{.info,.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}`.

Duplications and losses should be inferred using the locus tree and locus reconciliation.



## 2.3 Convert between reconciliation formats

### Convert from three-tree to LCT

    dlcpar convert \
        --3t_to_lct \
        -s config/paper.stree \
        -S config/paper.smap \
        data/paper/0/0.coal.tree

This creates the files `0.{lct.tree,lct.recon,lct.order}`.

### Convert from LCT to three-tree

    dlcpar convert \
        --lct_to_3t \
        -s config/paper.stree \
        -S config/paper.smap \
        -O .test \
        data/paper/0/0.lct.tree

This creates the files `0.test.{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters}`. Here, we use the output extension `-O .test` so as not to overwrite the original three-tree files.



## 2.4 Infer events from reconciliations

### Infer from LCT

    dlcpar events \
        --lct \
        -s config/paper.stree \
        -S config/paper.smap \
        data/paper/0/0.lct.tree

### Infer from three-tree

    dlcpar events \
        --3t \
        -s config/paper.stree \
        -S config/paper.smap \
        data/paper/0/0.coal.tree



## 2.5 Compare reconciliations

### Compare three-trees

Suppose we want to check for equality between two reconciliations in three-tree format. These reconciliations must map the same gene tree within the same species tree using the same species map.

    dlcpar equal \
        --3t \
        -s config/paper.stree \
        -S config/paper.smap \
        data/paper/0/0.coal.tree \
        data/paper/0/0.test.coal.tree

### Compare LCTs

Suppose we want to check for equality between two reconciliations in LCT format. These reconciliations must map the same gene tree within the same species tree using the same species map.

    dlcpar equal \
        --lct \
        -s config/paper.stree \
        -S config/paper.smap \
        data/paper/0/0.lct.tree \
        data/paper/0/0.dlcdp.lct.tree



## 2.6 Find MPR landscapes

### Find equivalent regions and events per region

    dlcpar landscape \
        -s config/paper.stree \
        -S config/paper.smap \
        -I .coal.tree -O .dlcscape \
        -x1234 --events I \
        data/paper/0/0.coal.tree

This creates the files `0.dlcscape{.info,.regions,.events}`.

### View landscape

    dlcpar view_landscape \
        data/paper/0/0.dlcscape.regions



## Clean up

To clean up all generated files in this example:

    find data/paper \
        -name '*dlc*' \
        -or -name '*lct* \
        -or -name '*test*' | \
        xargs rm

<!-- /Usage -->
