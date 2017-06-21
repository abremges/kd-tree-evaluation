# kd-tree-evaluation

## Download evaluation data

### ONT data
Source: http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/
```
wget https://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta
wget https://nanopore.s3.climb.ac.uk/MAP006-2_2D_pass.fasta
cat MAP006*.fasta | awk '/^>/{print ">" ++i; next}{print}' > ont.fasta
```

### PB data
Source: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly
```
wget https://s3.amazonaws.com/files.pacb.com/datasets/secondary-analysis/e-coli-k12-P6C4/p6c4_ecoli_RSII_DDR2_with_15kb_cut_E01_1.tar.gz
tar xfv p6c4_ecoli_RSII_DDR2_with_15kb_cut_E01_1.tar.gz
for f in E01_1/Analysis_Results/m141013*bax.h5; do pls2fasta -trimByRegion -minSubreadLength 1000 $f > $f.fasta; done
cat E01_1/Analysis_Results/m141013*.fasta | awk '/^>/{print ">" ++i; next}{print}' > pacbio.fasta
```

## Parse overlap files

Create a TSV with 2 fields: R1 and R2 (that overlap, excluding self-overlaps, ID(R1) < ID(R2)):  
``for f in *.out; do awk '{if ($1!=$2) print ($1<$2) ? $1"\t"$2 : $2"\t"$1}' < $f | sort | uniq > $f.tsv; done``
