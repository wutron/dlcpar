# 1. Usage

Running `dlcpar -h` will print out its command-line usage.

DLcpar has several sub-commands:  
| command | description |
|---|---|
| `convert`        | Convert between reconciliation structures|
| `equal`          | Check for equality between reconciliation structures |
| `events`         | Infer events in a reconciliation |
| `dp`             | Solve the MPR problem using dynamic programming |
| `search`         | Solve the MPR problem using heuristic search |
| `landscape`      | Find MPR landscapes across ranges of event counts |
| `view_recon`     | View LCT |
| `view_landscape` | View cost landscape |

# 2. File Formats

## 2.1 Reconciliations

Reconciliations are stored either in Labeled Coalescent Tree (LCT) or three-tree (3T) format. Most programs use the LCT, though you should use the 3T format if you need the locus tree in its own file.

### 2.1.1 Labeled Coalescent Tree (LCT)
| file | description |
|---|---|
| `X.lct.tree`  | coalescent tree (newick)
| `X.lct.recon` | species map and locus map |
| `X.lct.order` | partial order |

The species map and locus map file format is tab-delimited, where each line has three fields:
1. coal node ID
2. species node ID
3. locus

The partial order file format is tab-delimited, where each line has three fields:
1. species node ID
2. parent locus
3. list of coal node IDs (comma-separated)

### 2.1.2 Three-Tree

The three-tree format was originally developed for the [DLCoal](http://compbio.mit.edu/dlcoal) model. Further details can be found in the [documentation](http://compbio.mit.edu/dlcoal/pub/dlcoal/doc/dlcoal-manual.html) for the software package.

| file | description |
|---|---|
| `X.coal.tree`   | coalescent tree (newick) |
| `X.coal.recon`  | reconciliation between coal tree and locus tree |
| `X.locus.tree`  | locus tree (newick) |
| `X.locus.recon` | reconciliation between locus tree and species tree |                     
| `X.daughters`   | daughter nodes in locus tree (one per line) |

The reconciliation file format is tab-delimited, where each line has three fields:
1. "inner tree" node ID
2. "outer tree" node ID
3. event ("gene", "spec", "dup", "none")

This format allows for the coal tree to be mapped "inside" the locus tree or the locus tree "inside" the species tree. For the reconciliation between coal tree and locus tree, all events are "none". For the reconciliation between locus tree and species tree, events can be "gene", "spec", or "dup".

## 2.2 Landscapes

There exist two different formats for landscapes, one for regions and one for events.
| file | description |
|---|---|
| `X.regions` | representation of position and area of cost regions |
| `X.events`  | event count vector followed by counts of optimal reconciliations |

Landscape are stored as `X.regions`. Each line has a `shapely` region object with the position the region and the area of the region. Use `view_regions` to display the cost landscape for a region file. Each of the descriptions of the region begins with an event count vector represented in the form `<duplication, loss, coalescence>: # of solutions`. This is followed by the position of the region and then the area of the region. 

Events are stored as `X.events`. For the first part of the file, each line starts with an event count vector, followed by the count of optimal reconciliations with that count vector and the events in all of those reconciliations.

The event count vector is represented in the form  `<D,L,C>:#` representing the counts of duplicationes, losses, and coalescences, followed by the number of solutions.

Each event is indicated through a three-tuple, where the first element is the event type, the second element is a string representation of the event, and the third element is the species in which the event occured:


| event | symbol | string | species |
|---|---|---|---|
| `speciation` | S | (leaves on first subtree in alphabetical order )<br> (leaves on other subtree) | species above the speciation |
| `duplication`  | D |(leaves on subtree with daughter locus) <br>   (leaves on subtree with parent locus)| species where the dup occurred
| `loss`  | L | (leaves on other subtree that is not lost) |      species where the locus was lost
| `coal_spec` | C| (tuples of leaves for each subtree which did not coalesce) | species in which lineages could have coalesced           
| `coal_dup`   | K | (tuples of leaves for each subtree which did not coalesce at a dup) | species with the duplication

Each subtree is encased in parenthesis `()`.

For the second part of the file, each line starts with a number, followed by a list of events appearing in that number of regions. The list of events uses the same notation as above.



# 3. Programs

DLCpar uses the labeled coalescent tree (LCT) to infer the species and locus to which a gene belongs. The two main programs are DLCpar and DLCscape. The DLCpar program generates the most parsimonious reconciliation given the costs for the three event types, while the DLCscape program divides the cost space into regions. In each region, the total event counts for each reconciliation are the same.


## 3.1 Reconciliations

### 3.1.1 dp

The &nbsp; `dlcpar dp` &nbsp; program finds a most parsimonious gene tree-species tree reconciliation.

The simplest way to use the program:

    dlcpar dp -s <species tree> -S <species map> <gene tree>

where

    <gene tree>    - filename of gene tree
    <species tree> - filename of species tree
    <species map>  - filename of gene to species map

If the gene tree has filename &nbsp; `<X.tree>` , this returns a MPR in LCT format as &nbsp; `X.lct.tree` , &nbsp; `X.lct.recon` , &nbsp; `X.lct.order`.

#### Setting a seed or a log file

By default, the random number seed is randomly generated, and the program does not output a debugging log.   
To set a random number seed or to set to output a debugging log, add the following to the  command:

    -x, --seed <random seed>
    -l, --log  

where 

    <random seed> - random number seed

#### Changing file extensions

By default,&nbsp;  `dlcpar dp` &nbsp; has a blank input file extension, and output file extension is set to .dlcpar. To change the input/output file extensions, add the following to the command:

    -I, --inputext <input file extension>
    -O, --outputext <output file extension>

where 

    <input file extension>  - input file extension name
    <output file extension> - output file extension name

#### Changing output file format
By default,&nbsp;  `dlcpar dp`&nbsp;  outputs in LCT format. To output in 3-tree format, add the following to the command:

    --output_format {lct or 3t} 

If the command is followed by lct, the output format will be in lct.  
If the command is followed by 3t, the output format will be in 3-tree format.

#### Incorporating multiple samples per species
By default, the &nbsp; `dlcpar dp`&nbsp;  expects one sample per species. To encorporate a locus map, add the following to the command:

    --lmap <locus map>

where

     <locus map> - gene to locus map (species-specific)


#### Sampling multiple optimal reconciliations
By default, &nbsp; `dlcpar dp`&nbsp;  returns a single MPR sampled uniformly at random. To sample multiple optima (also uniformly at random), add the following parameter:

    -n <number of reconciliations>

This returns multiple MPRs in the same LCT format, but each file will contain the number of solutions as selected in the parameter, separated by '# Solution 0', '# Solution 1', etc. 

Note that `dlcpar search` &nbsp; cannot sample multiple optima. Furthermore, other programs (`dlcpar convert`, &nbsp; `dlcpar events`) expect a single solution per file, so you will have to separate the solutions into individual files to use those programs.

#### Changing event costs
By default, the reconciliation cost module used by &nbsp; `dlcpar dp`&nbsp;  assumes equal costs (D=1, L=1) for inferred (duplication-loss) events. To change this, add the following commands:

    -D, --dupcost <dup cost> 
    -L, --losscost <loss cost>   
    -C, --coalcost <coal cost>  
    -K, --coaldupcost <coal dup cost> 
where

    <dup cost>      - dupduplication cost (default: 1.0)
    <loss cost>     - loss cost (default: 1.0)
    <coal cost>     - deep coalescence cost at speciation (default: 0.5)
    <coal dup cost> - deep coalescence cost at duplication if different


#### Heuristics Settings:

The following prescreen settings are available:

    --no_prescreen
set to disable prescreen of locus maps


    --prescreen_min <prescreen min> 
prescreen locus maps if min (forward) cost exceeds this value (default: 50)

    --prescreen_factor <prescreen factor>  
prescreen locus maps if (forward) cost exceeds this factor x min (forward) cost (default: 10)

By defeault, the there is no maximum number of co-existing loci (per species), a maximum of 4 duplications (per species) and a maximum of 4 losses (per species). To set maximums, add the following command:   

    --max_loci <max # of loci> 
    --max_dups <max # of dups>   
    --max_losses <max # of losses>   

where

    <max # of loci>   - max # of co-existing loci (per species), set to -1 for no limit (default: -1)
    <max # of dups>   - max # of duplications (per species), set to -1 for no limit (default: 4)
    <max # of losses> - max # of losses (per species), set to -1 for no limit (default: 4)


### 3.1.2 search

 For some gene trees that are very large or highly incongruent to the species tree, you may need to use DLCpar-search, which relies on the 3T model and searches the space of reconciliations using a hill-climbing approach.

The simplest way to use the program:

    dlcpar search -s <species tree> -S <species map> <gene tree>

where

    <gene tree>    - filename of gene tree
    <species tree> - filename of species tree
    <species map>  - filename of gene to species map


If the gene tree has filename &nbsp; `<X.tree>`, this returns a MPR in 3T format as &nbsp; `X.coal.tree`,&nbsp;  `X.coal.recon`, &nbsp; `X.locus.tree`,&nbsp; `X.locus.recon`,&nbsp; `X.daughters`. However, this cannot be converted to LCT format because the locus tree is undated. Duplications and losses should be inferred using the locus tree and locus reconciliation.

#### Setting a seed or a log file

By default, the random number seed is randomly generated, and the program does not output a debugging log.   
To set a random number seed or to set to output a debugging log, add the following to the  command:

    -x, --seed <random seed>
    -l, --log  

where 

    <random seed> - random number seed

#### Changing file extensions

By default,&nbsp;  `dlcpar search` &nbsp; has a blank input file extension, and output file extension is set to .dlcpar. To change the input/output file extensions, add the following to the command:

    -I, --inputext <input file extension>
    -O, --outputext <output file extension>

where 

    <input file extension>  - input file extension name
    <output file extension> - output file extension name

#### Changing event costs
By default, the reconciliation cost module used by &nbsp; `dlcpar search`&nbsp;  assumes equal costs (D=1, L=1) for inferred (duplication-loss) events. To change this, add the following commands:

    -D, --dupcost <dup cost> 
    -L, --losscost <loss cost>   
    -C, --coalcost <coal cost>  

where

    <dup cost>  - dupduplication cost (default: 1.0)
    <loss cost> - loss cost (default: 1.0)
    <coal cost> - deep coalescence cost at speciation (default: 0.5)

#### search specific settings

By default, `dlcpar search` &nbsp; has the number of search iterations and prescreening iterations set to 10 and 20 respectively. This may be changed through adding the following commands:

    -i, --iter <# iterations>
    --nprescreen <# prescreens>

where

    <# iterations> - number of search iterations (default: 10)
    <# prescreens> - number of prescreening iterations (default: 20)
  
In addition,  `dlcpar search` &nbsp; has the option to stop search after a convergence (a solution has converged if it has not changed for the specifiednumber of iterations) and the specification of the initial locus tree for searching. This is done through adding the following command:

    --nconverge <# converge>
    --init-locus-tree <tree file>

where

    <# converge> - the number of solutions required to determine convergence
    <tree file> - initial locus tree filename for searching

### 3.1.3 landscape

The &nbsp; `dlcpar landscape` &nbsp; program finds MPR landscapes across ranges of event.

The simplest way to use the program:

    dlcpar landscape -s <species tree> -S <species map> <gene tree>

where

    <gene tree>    - filename of gene tree
    <species tree> - filename of species tree
    <species map>  - filename of gene to species map

If the gene tree has filename &nbsp; `<X.tree>` , this returns region and landscape files as &nbsp; `X.regions` , &nbsp; `X.events`.

#### Event specific settings

By default, if --eventsflag is not present, the program will not compute events at all, which is faster if you just want the cost space. Otherwise, --events reports (U)nion or (I)ntersection of events. This is altered by adding the following command:

    --events <output events>

where

    <output events> - U or I to set to output (I)ntersection or (U)nion of events (default: None)

The program also has a function which draws regions to screen. This can be activated through adding the following command:

    --draw_regions




#### Setting a seed or a log file

By default, the random number seed is randomly generated, and the program does not output a debugging log.   
To set a random number seed or to set to output a debugging log, add the following to the  command:

    -x, --seed <random seed>
    -l, --log  

where 

    <random seed> - random number seed

#### Changing file extensions

By default,&nbsp;  `dlcpar landscape` &nbsp; has a blank input file extension, and output file extension is set to .dlcscape. To change the input/output file extensions, add the following to the command:

    -I, --inputext <input file extension>
    -O, --outputext <output file extension>

where 

    <input file extension>  - input file extension name
    <output file extension> - output file extension name

#### Incorporating multiple samples per species
By default, the &nbsp; `dlcpar landscape`&nbsp;  expects one sample per species. To encorporate a locus map, add the following to the command:

    --lmap <locus map>

where

     <locus map> - gene to locus map (species-specific)


#### Setting the cost range

The cost ranges in the program are normalised to coalescence cost = 1. By defeault, the cost range for both duplication and loss costs are set from 0.2 to 5. This can be changed by adding the following to the command: 

    -D, --duprange <dup range>
    -L, --lossrange <loss range>

where

    <dup range> - duplication cost range (default: 0.2-5)
    <loss range> - loss cost range (default: 0.2-5)


#### Heuristics Settings
By defeault, the there is no maximum number of co-existing loci (per species), a maximum of 4 duplications (per species) and a maximum of 4 losses (per species). To set maximums, add the following command:   

    --max_loci <max # of loci> 
    --max_dups <max # of dups>   
    --max_losses <max # of losses>   

where

    <max # of loci>   - max # of co-existing loci (per species), set to -1 for no limit (default: -1)
    <max # of dups>   - max # of duplications (per species), set to -1 for no limit (default: 4)
    <max # of losses> - max # of losses (per species), set to -1 for no limit (default: 4)

## 3.2 Utilities

### 3.2.1 convert
This program converts between reconcilation structures of LCT and 3-tree. 

The simplest way to use the program:

For conversion from LCT to 3T format:

    dlcpar convert --lct_to_3tree -s <species tree> -S <species map>
    <lct tree> 

where

    <lct tree>     - filename of lct tree (needs to be with .lct.tree)
    <species tree> - filename of species tree
    <species map>  - filename of gene to species map

If the lct tree has filename &nbsp; `X.lct.tree` , this returns files of 3T format as &nbsp; `X.coal.tree` , &nbsp; `X.coal.recon` , &nbsp; `X.locus.tree`, &nbsp; `X.locus.recon`, &nbsp; `X.daughters`.

For conversion from 3T to LCT format:

    dlcpar convert --3tree_to_lct -s <species tree> -S <species map>
    <gene tree> 

where

    <gene tree>    - filename of gene tree (needs to be with .coal.tree)
    <species tree> - filename of species tree
    <species map>  - filename of gene to species map

If the gene tree has filename &nbsp; `X.coal.tree` , this returns files of LCT format as &nbsp; `X.lct.tree` , &nbsp; `X.lct.order` , &nbsp; `X.lct.recon`.

#### Miscellaneous settings (for 3tree_to_lct)

    -use-locus-recon
If set, use locus recon file rather than MPR files

    -no-delay 
If set, disallow duplication between speciation and coalescence.

### 3.2.2 equal

This program checks for equality of reconciliation structues. If true, it returns "True". Otherwise, it will result in a key error.

The simplest way to use the program:

    dlcpar equal -- format (3t or lct) -s <speices tree> -S <species map> 
     <prefix 1> <prefix 2>

where

    (3t or lct)    - the format to compare equality in
    <species tree> - filename of species tree
    <species map>  - filename of gene to species map
    <prefix 1>     - prefix of tree to be compared #1
    <prefix 2>     - prefix of tree to be compared #2

#### Miscellaneous settings (for 3t format)

    --use-locus-mpr      
 if set, use MPR rather than locus recon file

### 3.2.3 events 
Infer events in reconciliations.

The simplest way to use the program:

    dlcpar events --format (3t or lct) -s <speices tree> -S <species map>  -I, --inputext <input file extension> <gene tree>

where

    (3t or lct)            - the format to compare equality in
    <species tree>         - filename of species tree
    <species map>          - filename of gene to species map
    <input file extension> - input file extension (.coal.tree for 3t and .lct.tree for lct)
    <gene tree>            - filename of gene tree

#### Incorporating multiple samples per species
By default, the &nbsp; `dlcpar landscape`&nbsp;  expects one sample per species. To encorporate a locus map, add the following to the command:

    --lmap <locus map>

where

     <locus map> - gene to locus map (species-specific) (only for lct format)

#### Miscellaneous

    --by-fam (Outputs by family)
    --use-famid (Uses family id for the genes)
    --use-locus_recon  (if set, use locus recon rather than MPR [only for 3t format])

### 3.2.4 view_recon
Visualise a LCT.

The simpelest way to use the program:

    dlcpar view_recon -s <species tree> -o <output file> <prefix>

where

    <species tree> - filename of species tree
    <output file>  - filename of output file (.svg)
    <prefix> - prefix of the lct files to visualise

#### Visualisation settings
By default, either the output option (-o) or the viewer option (-v) must be selected. The following options may also be chosen:

    -v, --viewer <svg viewer>
        --xscale <x-scaling>
        --yscale <y-scaling>
        --names               display internal node names
        --snames              display species names

where

    <svg viewer> - name of the svg viewer
    <x-scaling>  - scaling in the x direction in the visualisation
    <y-scaling>  - scaling in the y direction in the visualisation


### 2.2.5 view_landscape
Visualises landscape.

The simplest way to use the program:

    dlcpar view_landscape <landscape>

where

    <landscape> - filename of the regions file created through landscape (.regions)

#### Visualisation settings
The program may also be set to output by adding the following to the command:

    -o <output>

where

    <output> - filename of the output (.pdf)

In additon, the logging may be seen through adding the following to the command:

    --linear



# 4. Examples

See `examples/test.sh` for an example of how to use each program in this package.
