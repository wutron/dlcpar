<h1>Table of Contents</h1>

<!-- TOC -->

- [1. Install](#1-install)
    - [1.1 Requirements](#11-requirements)
    - [1.2 Install](#12-install)
- [2. Overview](#2-overview)
    - [2.1 Main Programs](#21-main-programs)
    - [2.2 Utilities](#22-utilities)
    - [2.3 Visualizations](#23-visualizations)
- [3. File Formats](#3-file-formats)
    - [3.1 Trees and Species Maps](#31-trees-and-species-maps)
    - [3.2 Reconciliations](#32-reconciliations)
        - [Labeled Coalescent Tree (LCT)](#labeled-coalescent-tree-lct)
        - [Three-Tree](#three-tree)
    - [3.2 Landscapes](#32-landscapes)
        - [Regions](#regions)
        - [Events](#events)
- [4. Programs](#4-programs)
    - [4.1 Main Programs](#41-main-programs)
        - [4.1.1 dp](#411-dp)
            - [Setting the seed, log, or output file format](#setting-the-seed-log-or-output-file-format)
            - [Changing the file extensions](#changing-the-file-extensions)
            - [Incorporating multiple samples or sampling multiple solutions](#incorporating-multiple-samples-or-sampling-multiple-solutions)
            - [Changing the event costs](#changing-the-event-costs)
            - [Setting parameters for heuristic screening](#setting-parameters-for-heuristic-screening)
        - [4.1.2 ilp](#412-ilp)
            - [Setting the ILP solver and parameters](#setting-the-ilp-solver-and-parameters)
        - [4.1.3 search](#413-search)
            - [Setting search parameters](#setting-search-parameters)
        - [4.1.4 landscape](#414-landscape)
            - [Outputting events](#outputting-events)
            - [Visualizing regions](#visualizing-regions)
            - [Changing the event cost range](#changing-the-event-cost-range)
    - [4.2 Utilities](#42-utilities)
        - [4.2.1 convert](#421-convert)
            - [Changing the file extensions](#changing-the-file-extensions-1)
            - [Setting miscellaneous options](#setting-miscellaneous-options)
        - [4.2.2 equal](#422-equal)
            - [Changing the file extensions](#changing-the-file-extensions-2)
            - [Setting miscellanous options](#setting-miscellanous-options)
        - [4.2.3 events](#423-events)
            - [Setting miscellaneous options](#setting-miscellaneous-options-1)
    - [4.3 Visualizations](#43-visualizations)
        - [4.3.1 view_lct](#431-view_lct)
            - [Setting output options](#setting-output-options)
        - [4.3.2 view_landscape](#432-view_landscape)
            - [Setting output options](#setting-output-options-1)

<!-- /TOC -->





<!-- Install -->

<a id="install"></a>
# 1. Install



<a id="install-reqs"></a>
## 1.1 Requirements

This package has the following requirements:

- [Python](http://python.org/) (2.7.x)
- [Numpy](http://www.numpy.org/) (1.13 or greater)
- [Shapely](https://pypi.python.org/pypi/Shapely) (1.5 or greater) -- only for `landscape` and `view_landscape`
- [Matplotlib](http://matplotlib.org/) (1.5 or greater) -- only for `landscape` and `view_landscape`
- [PuLP](http://coin-or.github.io/pulp/) (2.1) -- only for `ilp`
- [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) (12.9) -- only for `ilp` (defaults to open-source CBC solver if not found)


<a id="install-install"></a>
## 1.2 Install

Run

```console
python setup.py install
```

If you do not have permission to install software on your system, you can install into another directory using the `--prefix` or `--home` flags to `setup.py`.

For example
```console
python setup.py install --prefix=/home/username/python
```

or
```console
python setup.py install --home=~
```

If you did not install in the standard bin directory, you will need to set your `PATH` variable to the alternate location.

If you did not install in the standard Python site-packages directory, you will need to set your `PYTHONPATH` variable to the alternate location. See [Python documentation](https://docs.python.org/2/using/cmdline.html#environment-variables) for further details.

<!-- /Install -->





<!-- Overview -->

<a id="overview"></a>
# 2. Overview

Running `dlcpar -h` will print out its command-line usage.

DLCpar has several commands used in various situations.

## 2.1 Main Programs
| command | description |
|---|---|
| [`dp`](#progs-main-dp)                         | Solve the MPR problem using dynamic programming |
| [`ilp`](#progs-main-ilp)                       | Solve the MPR problem using integer linear programming |
| [`search`](#progs-main-search)                 | Solve the MPR problem using heuristic search |
| [`landscape`](#progs-main-landscape)           | Find MPR landscapes across range of event costs |

## 2.2 Utilities
| command | description |
|---|---|
| [`convert`](#progs-utils-convert)               | Convert between reconciliation structures |
| [`equal`](#progs-utils-equal)                   | Check for equality between reconciliation structures |
| [`events`](#progs-utils-events)                 | Infer event counts in a reconciliation |

## 2.3 Visualizations
| command | description |
|---|---|
| [`view_lct`](#progs-vis-viewlct)             | View the labeled coalescent tree |
| [`view_landscape`](#progs-vis-viewlandscape) | View landscape of equivalent regions |

See [examples](examples/EXAMPLES.md) of how to use each command in this package.

<!-- /Overview -->





<!-- File Formats -->

<a id="formats"></a>
# 3. File Formats



<a id="formats-trees"></a>
## 3.1 Trees and Species Maps

Trees (`*.tree`, `*.stree`) should be specified using the Newick file format.

There are several restrictions on what IDs are allowed. Many of these restrictions are common for other similar phylogenetic software. The safest IDs follow these restrictions:

1. the first and last characters of the ID are a-z A-Z 0-9 _ -
2. the middle characters can be a-z A-Z 0-9 _ - . or the space character ' '
3. the ID should not be purely numerical characters 0-9
4. the ID should be unique within each tree

Space characters are discouraged since they will probably cause problems with other bioinfomatics software that you may use. Characters such as parentheses "(" ")" and colons ":" are not allowed because the newick file format uses these characters for describing the structure of the tree.

It is also easier to use gene IDs that have a prefix or suffix that indicates the species ID. For example "human_HOXC5" is a human gene. This is not a requirement, but it does make preparing a gene to species mapping file (`*.smap`) easier.

If IDS are not given to the ancestral nodes, the program will by default name them with "nXXX" where XXX is determined via preorder traversal of the tree.


Species maps (`*.smap`) should be specified in a tab-delimited format, where each line has two fields:

1. pattern matching a gene ID
2. species ID

Only 3 types of gene ID patterns are supported. The pattern can either be an exact matching string, a prefix (denoted "text\*"), or a suffix (denoted "\*text"). The "*" is the only special wildcard character.

The species ID should be the same as those used in the species tree. All patterns and IDs are case-sensitive.

When mapping a gene ID to a species ID, all exact matches are processed first. If no exact match is found, the patterns are then processed in the same order as they appear in the file until a match is found.



<a id="formats-recons"></a>
## 3.2 Reconciliations

Reconciliations are represented either in Labeled Coalescent Tree (LCT) or three-tree (3T) format. Most programs use the LCT, though you should use the 3T format if you need the locus tree in its own file. Unless otherwise specified, the lines can be given in any order within the file.

<a id="formats-recons-lct"></a>
### Labeled Coalescent Tree (LCT)

The reconciliation is represented as a tuple of three variables: the species map *M*, the locus map *L*, and the partial order *O*.  Each of these variables is stored in one of the following output files.

| file | description |
|---|---|
| `*.lct.tree`  | coalescent tree |
| `*.lct.recon` | species map and locus map |
| `*.lct.order` | partial order |

The coalescent tree is specified in Newick file format.  It is a copy of the gene tree with internal names and not technically part of the reconciliation. This file is created in case the original input gene tree did not specify internal node names.

The species map and locus map is tab-delimited, where each line has three fields:
1. coal node ID
2. species node ID
3. locus

Each line specifies the mapping of one node in the gene tree (field 1) to one node in the species tree (field 2) and to one locus (field 3).

The partial order is tab-delimited, where each line has three fields:
1. species node ID
2. parent locus
3. list of coal node IDs (comma-separated)

For a species node (field 1) and locus (field 2) that duplicates, each line specifies the order of nodes (field 3) that map to that species and locus or that map to that species and descend from that locus.

Use [`view_lct`](#progs-vis-viewlct) to view the reconciliation.

<a id="formats-recons-3t"></a>
### Three-Tree

The reconciliation is represented as a tuple of four variables: the locus tree *T<sup>L</sup>*, the gene-to-locus tree mapping *R<sup>G</sup>*, the locus-to-species mapping *R<sup>L</sup>*, and the set of daughter nodes *&delta;<sup>L</sup>*. Each of these variables is stored in one of the following output files.

| file | description |
|---|---|
| `*.coal.tree`   | coalescent tree (newick) |
| `*.coal.recon`  | reconciliation between coal tree and locus tree |
| `*.locus.tree`  | locus tree (newick) |
| `*.locus.recon` | reconciliation between locus tree and species tree |
| `*.daughters`   | daughter nodes in locus tree (one per line) |

The coalescent tree is specified in Newick file format. It is a copy of the gene tree with internal names and not technically part of the reconciliation. This file is created in case the original input gene tree file did not specify internal node names.

The reconciliation file format is tab-delimited, where each line has three fields:
1. "inner tree" node ID
2. "outer tree" node ID
3. event ("gene", "spec", "dup", "none")

This format allows for the coalescent tree to be mapped "inside" the locus tree or the locus tree "inside" the species tree. Each line specifies the mapping of one node in the inner tree (field 1) to one node or branch in the outer tree (field 2). Branches are indicated using the node ID directly below it (i.e. the younger of the two incident nodes). For the reconciliation between coalescent tree and locus tree, all events are "none". For the reconciliation between locus tree and species tree, events can be "gene", "spec", or "dup".

The daughters file format is simple, with each line containing one locus node ID that belongs to the set of daughter nodes in the locus tree.

The three-tree format was originally developed for the [DLCoal](http://compbio.mit.edu/dlcoal) model. More details can be found in the [documentation](http://compbio.mit.edu/dlcoal/pub/dlcoal/doc/dlcoal-manual.html#sec-file-recon) for the software package.



<a id="formats-landscapes"></a>
## 3.2 Landscapes

Landscapes explore Pareto-optimal reconciliations and events within these reconciliations.

| file | description |
|---|---|
| `*.regions` | equivalent regions |
| `*.events`  | events per region |

<a id="formats-landscapes-regions"></a>
### Regions

The regions file format is tab-delimited. The first line specifies the input cost ranges *d<sub>min</sub>*, *d<sub>max</sub>*, *l<sub>min</sub>*, *l<sub>max</sub>*. Each subsequent line specifies one region and has three fields:
1. event count descriptor
2. region type
3. coordinates
4. area

Event count descriptors (field 1) are in the form `<d,l,c>:#`, where `d` is the number of duplications, `l` the number of losses, `c` the number of coalescences, and `#` the number of reconciliations with event count `<d,l,c>`. The region type (field 2) is either `POLYGON`, `LINE`, or `POINT` corresponding to a [`shapely`](https://pypi.python.org/pypi/Shapely) region object. The coordinates (field 3) and area (field 4) specify the coordinates and area of the region.

Use [`view_landscapes`](#progs-vis-viewlandscape) to view the landscape.

<a id="formats-landscapes-events"></a>
### Events

The events file format is tab-delimited and has two parts.

In the first part of the file, each line specifies one region and has several fields. Field 1 specifies the event count descriptor (see [above](#formats-landscapes-regions)) for the region. Subsequent fields specify the set of all events across all reconciliations with that event count.

An event is represented as a tuple, where the first element is the event type, the second element is a string representation of the event, and the third element is the species in which the event occurred. Because internal node names can be arbitrary, events are represented using subtrees of the gene tree (for coalescence events) or subtrees of the locus tree (for speciation, duplication, and loss events). String representations are based on the gene relationship file of DLCoal (see [documentation](http://compbio.mit.edu/dlcoal/pub/dlcoal/doc/dlcoal-manual.html#sec-file-rel)).

| event | symbol | description | species |
|---|---|---|---|
| `speciation`  | S | leaves on first subtree in alphabetical order<br> leaves on other subtree     | species above the speciation |
| `duplication` | D | leaves on subtree with daughter locus<br> leaves on subtree with parent locus | species where the duplication occurred |
| `loss`        | L | leaves on subtree that is not lost                                            | species where the locus was lost |
| `coal_spec`   | C | tuples of leaves for each subtree which did not coalesce                      | species in which lineages could have coalesced |
| `coal_dup`    | K | tuples of leaves for each subtree which did not coalesce at a duplication     | species with the duplication |

Each subtree is encased in parenthesis `()`. As with event count descriptors, an event may be followed by the number of reconciliations with that event, i.e. `string:#`.

In the second part of the file, each line summarizes the events that occur in exactly *k*, exactly *k-1*, ..., exactly *1* region. Field 1 specifies a number of regions. Subsequent fields specify the set of all events that occur in exactly that number of regions.

<!-- /File Formats -->





<!-- Programs -->

<a id="progs"></a>
# 4. Programs



<a id="progs-main"></a>
## 4.1 Main Programs

This software package has several main goals:
1. Infer one or more MPRs given a single setting of event costs ([`dp`](#progs-main-dp) and [`search`](#progs-main-search))
2. Partition the space of event costs into regions such that all event costs within the same region yield the same set of MPRs ([`landscape`](#progs-main-landscape))

Each main program will output an information file (`*.info`) that logs the run. This file includes the program version, executed command, random number generator seed, and runtime. There may also be other information depending on the program.


<a id="progs-main-dp"></a>
### 4.1.1 dp

The `dp` command finds a most parsimonious gene tree-species tree reconciliation through dynamic programming.

The simplest way to use the program is as follows:
```console
dlcpar dp -s <species tree> -S <species map> <gene tree>
```
where
```console
<species tree> - filename of species tree
<species map>  - filename of gene to species map
<gene tree>    - filename of gene tree
```

If the gene tree has filename `X`, this program returns a reconciliation in [LCT format](#formats-recons-lct) as `X.dlcdp.lct.tree`, `X.dlcdp.lct.recon`, `X.dlcdp.lct.order`.

#### Setting the seed, log, or output file format

```console
-x, --seed <random seed>
```
This sets a seed for the random number generator, allowing one to produce deterministic results. If not specified, the seed of the random number generator is based on the clock time.

```console
-l, --log
```
This outputs a debugging log to `X.dlcdp.log.gz`.

```console
--output_format {(lct)|3t}
```
This specifies the [output format](#formats-recons) for the reconciliation (default=`lct`). If `format=3t`, then a gene tree with filename `X` returns a reconciliation in [three-tree format](#formats-recons-3t) as `X.dlcdp.coal.tree`, `X.dlcdp.coal.recon`, `X.dlcdp.locus.tree`, `X.dlcdp.locus.recon`, `X.dlcdp.daughters`.

#### Changing the file extensions

```console
-I, --inputext <input file extension>
-O, --outputext <output file extension>
```
This specifies the input file extension (default=`''`) of the gene tree and the output file extension (default=`.dlcdp`) for all output files. The input extension will be stripped from the input filename to determine the output file prefix `PREFIX`. Every output file will have the form `$PREFIX$OUTPUTEXT.*` (e.g. `prefix.dlcdp.lct.tree`).

#### Incorporating multiple samples or sampling multiple solutions

By default, the program expects one sample per species and samples a single MPR uniformly at random.

Multiple samples are allowed by specifying a locus map file.
```console
--lmap <locus map>
```
The file format is tab-delimited, where each line has two fields:
1. pattern matching a gene ID
2. species-specific locus

The locus map is similar to the [species map](#formats-trees) except that genes are mapped to loci rather than to species.

```console
-n <number of reconciliations>
```
This specifies the number of MPRs to sample uniformly at random (default=`1`). The reconciliation files will contain multiple solutions, each separated by "# Solution 0", "# Solution 1", etc. Note that all utility and visualization commands expect a single solution per file, so the multiple solutions will have to be separated to use these commands.

#### Changing the event costs

```console
-D, --dupcost <dup cost>
-L, --losscost <loss cost>
-C, --coalcost <coal cost>
-K, --coaldupcost <coal dup cost>
```
This specifies the costs for each type of event (defaults *D=1*, *L=1*, *C=0.5*, *K=C*).

#### Setting parameters for heuristic screening

By default, the method uses heuristics to prescreen locus maps. This limits the search space, which might be required for large gene trees or species trees. However, it may also cause the method to find a sub-optimal reconciliation.

```console
--no_prescreen
```
This disables prescreening of locus maps.

If prescreened is enabled, several settings can be specified.

```console
--prescreen_min <prescreen min>
--prescreen_factor <prescreen factor>
```
Locus maps are enumerated by traversing the species tree in pre-order. Maps are filtered if the minimum cost exceeds the minimum value (default=`50`) or if the cost exceeds the factor value (default=`10`) multipled by the minimum cost.

```console
--max_loci <max # of loci>
--max_dups <max # of dups>
--max_losses <max # of losses>
```
This specifies the maximum number of loci (default=`inf`), duplications (default=`4`), or losses (default=`4`) allowed per species.


<a id="progs-main-ilp"></a>
### 4.1.2 ilp

The `ilp` command finds a most parsimonious gene tree-species tree reconciliation through integer linear programming. It is a useful alternative to `dlcpar ilp` when (1) the gene tree is too large or is highly incongruent to the species tree, as these cases may be too complex for the dynamic programming approach, or (2) you wish to limit the maximum runtime or memory when inferring an MPR.

The command works similarly to [`dp`](#progs-main-dp) command though with fewer options. The simplest way to use the program is as follows:

```console
dlcpar ilp -s <species tree> -S <species map> <gene tree>
```
where
```
<species tree> - filename of species tree
<species map>  - filename of gene to species map
<gene tree>    - filename of gene tree
```

If the gene tree has filename `X`, this program returns a reconciliation in [LCT format](#formats-recons-lct) as `X.dlcilp.lct.tree`, `X.dlcilp.lct.recon`, `X.dlcilp.lct.order`.

#### Setting the ILP solver and parameters

```console
--solver {(CBC_CMD)|CPLEX_PY}
-t, --time_limit <time limit>
-m, --mem_limit <mem limit>
-T, --threads <number of threads>
```
This specifies the ILP solver to use (default=`CBC_CMD`), the time limit in seconds (default: no limit), the memory limit in MB (default: no limit), and the number of threads (default: solver default).  You should limit the number of threads to the number of available processors.


<a id="progs-main-search"></a>
### 4.1.3 search

The `search` command finds a most parsimonious gene tree-species tree reconciliation through a heuristic search. This method searches the space of reconciliations through hill-climbing (or, rather, hill-descending) similar to the [DLCoalRecon](http://compbio.mit.edu/dlcoal/) probabilistic reconciliation method. It is a "last-resort" alternative to `dlcpar dp` and `dlcpar ilp`.

The command works similarly to [`dp`](#progs-main-dp) command though with fewer options. The simplest way to use the program is as follows:

```console
dlcpar search -s <species tree> -S <species map> <gene tree>
```
where
```
<species tree> - filename of species tree
<species map>  - filename of gene to species map
<gene tree>    - filename of gene tree
```

If the gene tree has filename `X`, this program returns a reconciliation in [three-tree format](#formats-recons-3t) as `X.dlcsearch.coal.tree`, `X.dlcsearch.coal.recon`, `X.dlcsearch.locus.tree`, `X.dlcsearch.locus.recon`, `X.dlcsearch.daughters`.

See [`dp`](#progs-main-dp) for a description of shared arguments. The default output file extension is `.dlcsearch`.

#### Setting search parameters

```console
-i, --iter <# iterations>
--nprescreen <# prescreens>
```
This specifies the number of search iterations to perform (default=`10`) and the number of prescreening iterations to use in the locus tree search (default=`20`).

```console
--nconverge <# converge>
```
Setting a threshold for convergence will stop the search if the optimal solution has not changed for the specified number of iterations.

```console
--init-locus-tree <tree file>
```
This specifies the initial locus tree for the search. By default, the gene tree topology is used as the initial locus tree.


<a id="progs-main-landscape"></a>
### 4.1.4 landscape

The `landscape` command finds MPR landscapes across a range of event costs.

The simplest way to use the program is as follows:

```console
dlcpar landscape -s <species tree> -S <species map> <gene tree>
```
where
```console
<species tree> - filename of species tree
<species map>  - filename of gene to species map
<gene tree>    - filename of gene tree
```

If the gene tree has filename `X` , this program returns [regions](#formats-landscapes-regions) and [events](#formats-landscapes-events) files as `X.dlcscape.regions` and `X.dlcscape.events`.

See [`dp`](#progs-main-dp) for a description of shared arguments. The default output file extension is `.dlcscape`.

#### Outputting events

```console
--events <output events>
```
By default, the program outputs only regions and not events. Events are not computed, so the program executes faster. Events can be reported using either union (`U`) or intersection (`I`). For an event count, union tracks events found in *any* MPR with the event count, and intersection tracks events found in *all* reconciliations with that event count.

#### Visualizing regions

```console
--draw_regions
```
This specifies to draw regions to screen. The regions can also be visualized separately using [`view_landscape`](#progs-vis-viewlandscape).

#### Changing the event cost range

```console
-D, --duprange <dup low>:<dup high>
-L, --lossrange <loss low>:<loss high>
```
This specifies the range of costs for duplications (default=`0.2:5`) and losses (default=`0.2:5`), each relative to the unit cost of coalescence.



<a id="progs-utils"></a>
## 4.2 Utilities


<a id="progs-utils-convert"></a>
### 4.2.1 convert

The `convert` command converts reconciliations between the [LCT format](#formats-recons-lct) and the [three-tree format](#formats-recons-3t).

The simplest way to use the program is as follows:
```console
dlcpar convert <conversion> -s <species tree> -S <species map> <gene tree>
```
where
```console
<conversion>   - either --3t_to_lct or --lct_to_3t
<species tree> - filename of species tree
<species map>  - filename of gene to species map
<gene tree>    - filename of gene tree
```

If `--3t_to_lct` is specified, the gene tree must have extension `.coal.tree`. This program expects a reconciliation in three-tree format with filenames `X.coal.tree`, `X.coal.recon`, `X.locus.tree`, `X.locus.recon`, `X.daughters` and returns a reconciliation in LCT format as `X.lct.tree`, `X.lct.recon`, `X.lct.order`.

If `--lct_to_3t` is specified, the gene tree must have extension `.lct.tree`. This program expects a reconciliation in LCT format with filenames `X.lct.tree`, `X.lct.recon`, `X.lct.order` and returns a reconciliation in three-tree format as `X.coal.tree`, `X.coal.recon`, `X.locus.tree`, `X.locus.recon`, `X.daughters`.

#### Changing the file extensions

```console
-I, --inputext <input file extension>
-O, --outputext <output file extension>
```

This specifies the input file extension (default=`.coal.tree` for `--3t_to_lct`, default=`.lct.tree` for `--lct_to_3t`) of the gene tree and the output file extension (default=`''`) for all output files. The input extension will be stripped from the input filename to determine the prefix `PREFIX`. Every input file in the reconciliation must start with `$PREFIX`, and every output file will have the form `$PREFIX$OUTPUTEXT.*`.

#### Setting miscellaneous options

The following options are only used for `--3t_to_lct`.

```console
--use-locus-recon
```
By default, the program uses the lowest common ancestor (LCA) mapping between the locus tree and the species tree. This specifies to use the actual `*.locus.recon` file instead.

```console
--no-delay
```
By default, the program includes delayed speciation nodes. However, such nodes are not needed if no duplications occur between speciation and coalescence. This specifies to remove such nodes. The program throws an exception if a node cannot be removed.


<a id="progs-utils-equal"></a>
### 4.2.2 equal

The `equal` command checks for equality of reconciliation structures.

The simplest way to use the program is as follows:
```console
dlcpar equal <format> -s <species tree> -S <species map> <gene tree 1> <gene tree 2>
```
where
```console
<format>       - either --3t or --lct
<species tree> - filename of species tree
<species map>  - filename of gene to species map
<tree 1>       - filename of first gene tree
<tree 2>       - filename of second gene tree
```

If `--3t` is specified, the gene trees must have extension `.coal.tree`. This program expects two reconciliations in three-tree format, each with filenames `X.coal.tree`, `X.coal.recon`, `X.locus.tree`, `X.locus.recon`, `X.daughters`.

If `--lct` is specified, the gene tree must have extension `.lct.tree`. This program expects two reconciliations in LCT format, each with filenames `X.lct.tree`, `X.lct.recon`, `X.lct.order`.

The program prints `True` or `False` to standard output. If the two reconciliations are not equal, it also prints the component that mismatched to standard error.

#### Changing the file extensions

```console
-I, --inputext <input file extension>
```

This specifies the input file extension (default=`.coal.tree` for `--3t`, default=`.lct.tree` for `--lct`) of the gene tree. The input extension will be stripped from the input filename to determine the prefix `PREFIX`. Every input file in the reconciliation must start with `$PREFIX`, and every output file will have the form `$PREFIX$OUTPUTEXT.*`.

#### Setting miscellanous options

The following options are only usef for `--3t`.

```console
--use-locus-lca
```
By default, the program uses the actual `*.locus.recon` file for the mapping between the locus tree and the species tree. This specifies to use the LCA mapping instead.


<a id="progs-utils-events"></a>
### 4.2.3 events

The `events` command counts events in reconciliations.

The simplest way to use the program is as follows:
```console
dlcpar events <format> -s <species tree> -S <species map>  <gene tree>
```
where
```console
<format>               - either --3t or --lct
<species tree>         - filename of species tree
<species map>          - filename of gene to species map
<gene tree>            - filename of gene tree
```

This prints a tab-delimited table, where each line has the following fields:
1. species node ID
2. parent species node ID
3. species branch length
4. number of genes in the species branch
5. number of duplications in the species branch
6. number of losses in the species branch
7. number of extra lineages due to deep coalescences in the species branch
8. number of gene births in the species branch

#### Setting miscellaneous options

By default, the program aggregates counts across all gene families.

```console
--by-fam
--use-famid
```
`--by-fam` specifies to output counts individually for each gene family. The table includes a new (left-most) column indicating the family id. This id is filename of the gene tree by default or the name of bottom-most directory if `--use-famid` is set, e.g. `PATH/FAMID/BASENAME`.

```console
--use-locus-recon
```
By default, the program uses the lowest common ancestor (LCA) mapping between the locus tree and the species tree. This specifies to use the actual `*.locus.recon` file instead.



## 4.3 Visualizations


<a id="progs-vis-viewlct"></a>
### 4.3.1 view_lct

The `view_lct` command visualizes a LCT.

The simplest way to use the program is as follows:
```console
dlcpar view_lct -s <species tree> <gene tree>
```
where
```console
<species tree> - filename of species tree
<gene tree>    - filename of gene tree
```

The gene tree must have extension `.lct.tree`. This program expects a reconciliation in LCT format with filenames `X.lct.tree`, `X.lct.recon`, `X.lct.order`.

#### Setting output options

By default, the program draws the LCT in SVG format then displays it to screen.

```console
-v, --viewer <svg viewer>
```
This specifies the svg viewer (default=`display`).

```console
-o <output file>
```
This specifies the output file. The image is not displayed to screen.

The following options may also be specified:
```console
--xscale <x-scaling>
--yscale <y-scaling>
--names
--snames
```
This sets the x-scale factor (default=`50`) and the y-scale factor (default=`50`), and specifies whether to display internal node names and species names.


<a id="progs-vis-viewlandscape"></a>
### 4.3.2 view_landscape

The `view_landscape` command visualizes the landscape of equivalent regions.

The simplest way to use the program is as follows:
```console
dlcpar view_landscape <regions>
```
where
```console
<regions> - filename of regions file created through landscape command
```

The terminal is blocked until the figure is closed.

#### Setting output options

```console
-o <output>
```
This specifies the output file. The file format is inferred from the extension of the filename. Supported formats depend on the active `matplotlib` backend. Most backends support png, pdf, ps, eps and svg.

```
--linear
```
This specifies the use of linear axes (rather than the default log axes).

<!-- /Programs -->
