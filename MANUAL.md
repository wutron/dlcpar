# Usage

Running `dlcpar -h` will print out its command-line usage.

DLcpar has several sub-commands:  
| command | description |
|---|---|
| `convert` | Convert between reconciliation structures|
| `equal`   | Check for equality between reconciliation structures |
| `events`  | Infer events in a reconciliation |
| `dp`      | Solve the MPR problem using dynamic programming |
| `search`  | Solve the MPR problem using heuristic search |



# Programs

DLCpar uses the labeled coalescent tree (LCT) to infer the species and locus to which a gene belongs. The two main programs are DLCpar and DLCscape. The DLCpar program generates the most parsimonious reconciliation given the costs for the three event types, while the DLCscape program divides the cost space into regions. In each region, the total event counts for each reconciliation are the same.


## Reconciliations

### dp

### search

## Utilities

### convert

### equal

### events 



# File Formats

## Reconciliations

Reconciliations are stored either in Labeled Coalescent Tree (LCT) or three-tree (3T) format. Most programs use the LCT, though you should use the 3T format if you need the locus tree in its own file.

### Labeled Coalescent Tree (LCT)
| file | description |
|---|---|
| `X.tree`  | coalescent tree (newick)
| `X.recon` | species map and locus map |
| `X.order` | partial order |

The species map and locus map file format is tab-delimited, where each line has three fields:
1. coal node ID
2. species node ID
3. locus

The partial order file format is tab-delimited, where each line has three fields:
1. species node ID
2. parent locus
3. list of coal node IDs (comma-separated)

### Three-Tree

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

## Events


## Landscapes

Landscape are stored as `X.regions`. Each line has a `shapely` region object with the position the region and the area of the region. Use `view_regions` to display the cost landscape for a region file.

## Events over landscapes

Events are stored as `X.events`. For the first part of the file, each line starts with an event count vector, followed by the count of optimal reconciliations with that count vector and the events in all of those reconciliations. For the second part of the file, each line starts with a number, followed by a list of events appearing in that number of regions.



# Examples

See `examples/test.sh` for an example of how to use each program in this package.
