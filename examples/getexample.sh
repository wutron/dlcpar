#====================
# single sample per species

path=~/work/dlcpar

# get example from dlcpar project (sim-flies, 10e6, 1x, famid 0)
mkdir config
cp -L $path/config/flies{.stree,.smap} config/

mkdir -p sim-flies/0
cp -L $path/data/sim-flies/10e6-1x/0/0{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters} data/sim-flies/0/

#====================
# multiple samples per species

path=~/work/dlcpar2

# get example from dlcpar project (sim-flies, 10e6, 1x, famid 0)
cp -L $path/config/flies.lmap config/

cp -L $path/data/sim-flies/10e6-1x/0/0.sample2{.coal.tree,.coal.recon,.locus.tree,.locus.recon,.daughters} data/sim-flies/0/

