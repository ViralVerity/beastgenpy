# beastgenpy

A python package to make XMLs for BEASTX that are not possible and/or annoying to make in BEAUti. Alternative to beastgen written in python.

Specifically it's useful for making single XMLs with multiple tree models, although it will also work with a single tree model.


Currently implemented models:

Phylogeography:
- Continuous with or without providing polygons to deal with uncertain locations
- Discrete with or without GLM (coming soon)

Clock models:
- Standard UCLD
- Strict

Population models:
- HMC Skygrid (currently broken)
- Constant population

Substitution models:
- GTR with rate variation

Tree stuff:
- Estimating topology entirely
- Starting tree
- Empirical trees
- Fixed tree
- Tip date sampling (will happen automatically for sequences with only a year)
- Codon partitioning (needs fixing)

