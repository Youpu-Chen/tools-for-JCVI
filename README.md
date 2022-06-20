# tools for JCVI 

This repo is used to tentatively deposit some scripts I used in my synteny analysis using JCVI

> why "tentatively"? 
>
> Because I might continue adding new scripts or revising the old scripts.



### (1) within- or between species synteny analysis 

1. we need to create a configure file.

```shell
vim jcvi.config
# species1 species1
# species1 species2
```

As shown above and in this repo, the configure file contains the species you want to explore the synteny relationship.

Each line represents a analysis that conduct within species or between two species, 

e.g. 

- `species1 species1` means run JCVI on its own genome.
- `species1 species2`means run JCVI between species1 and species2.

2. After specified the configure file, now it's time to run the JCVI automatically, not manually.

```shell
sh run_jcvi.sh
```

Finally, there are some matters needing attention: The config file you specified in the `run_jciv.sh`, it should name after the config you created.



### (2) JCVI painting

`v2.jcvi_painting.py` is used to paint the synteny relationships without manually changing the `.simple` file which is created by the JCVI.

1. color config file

for example, there are 12 chromosomes/blocks in the genome you studied and you want to set 12 different colors for the each chromosomes/blocks, you should create a `color config` (this file is tab-delimited) like below:

```shell
1	#fa8072
2	#003153
3	#ffff4d
4	#800080
5	#6495ED
6	#00FFFF
7	#98FB98
8	#7FFF00
9	#808000
10	#FFFACD
11	#FFA500
12	#CD5C5C
```

Then, run the python scripts above,

```shell
python3 v2.jcvi_painting.py --bed species1.bed --config color.config --simple species1.species2.anchors.simple --output new.species1_species2.simple

# save the origin .simple file
mv species1.species2.anchors.simple species1.species2.anchors.simple.back
mv new.species1_species2.simple species1.species2.anchors.simple
# run JCVI macro visualization
python -m jcvi.graphics.karyotype seqids layout
```

