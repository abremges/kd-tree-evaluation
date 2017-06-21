# kd-tree-evaluation

## Download and prepare evaluation data

ONT data: http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/
```
wget https://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta
wget https://nanopore.s3.climb.ac.uk/MAP006-2_2D_pass.fasta
cat MAP006*.fasta | awk '/^>/{print ">" ++i; next}{print}' > ont.fasta
```

PB data: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly
```
wget https://s3.amazonaws.com/files.pacb.com/datasets/secondary-analysis/e-coli-k12-P6C4/p6c4_ecoli_RSII_DDR2_with_15kb_cut_E01_1.tar.gz
tar xfv p6c4_ecoli_RSII_DDR2_with_15kb_cut_E01_1.tar.gz
for f in E01_1/Analysis_Results/m141013*bax.h5; do pls2fasta -trimByRegion -minSubreadLength 1000 $f > $f.fasta; done
cat E01_1/Analysis_Results/m141013*.fasta | awk '/^>/{print ">" ++i; next}{print}' > pacbio.fasta
```

## Generate reference-based ground truth
Reference: https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.2 (safe as ``reference.fasta``)
```
bwa index reference.fasta
bwa mem -x ont2d reference.fasta ont.fasta | samtools view -F 2048 -uT reference.fasta - | samtools sort - > ont.fasta.bam
bwa mem -x pacbio reference.fasta pacbio.fasta | samtools view -F 2048  -uT reference.fasta - | samtools sort - > pacbio.fasta.bam
for f in *.bam; do bedtools bamtobed -i $f > $f.bed; done
for f in *.bed; do bedtools intersect -wo -a $f -b $f > $f.self; done
for f in *.self; do cut -f4,10 $f | awk '{if ($1!=$2) print ($1<$2) ? $1"\t"$2 : $2"\t"$1}' | sort | uniq > $f.tsv; done
```
``ont.fasta.bam.bed.self.tsv`` and ``pacbio.fasta.bam.bed.self.tsv`` are tab-separated files with only 2 fields: ID(R1) and ID(R2); R1 and R2 are two reads that overlap, ID(R1) < ID(R2). 

## Generate (and parse) ``kd`` read overlaps
Install ``kd``: https://github.com/dzif/kd-tree-overlapper
```
kd -o ont.kd.out -i ont.fasta
kd -o pacbio.kd.out -i pacbio.fasta
for f in *.out; do awk '{if ($1!=$2) print ($1<$2) ? $1"\t"$2 : $2"\t"$1}' < $f | sort | uniq > $f.tsv; done
```
Again, ``ont.kd.out.tsv`` and ``pacbio.kd.out.tsv`` are tab-separated files with only 2 fields: ID(R1) and ID(R2); R1 and R2 are two reads that overlap, ID(R1) < ID(R2). 

## Compare overlaps and generate result table
This could (and probably should) be scripted, but we'll do it manually:

