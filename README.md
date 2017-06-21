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

## Generate (and parse) read overlaps
Install ``kd``: https://github.com/dzif/kd-tree-overlapper
```
kd -o ont.kd.out -i ont.fasta
kd -o pacbio.kd.out -i pacbio.fasta
for f in *.out; do awk '{if ($1!=$2) print ($1<$2) ? $1"\t"$2 : $2"\t"$1}' < $f | sort | uniq > $f.tsv; done
```
Again, ``ont.kd.out.tsv`` and ``pacbio.kd.out.tsv`` are tab-separated files with only 2 fields: ID(R1) and ID(R2); R1 and R2 are two reads that overlap, ID(R1) < ID(R2). 

## Compute sensitivity, precision, and F1 score
This could – and probably should – be scripted, but we'll do it manually (default MHAP values in parentheses):

### ONT
P = `wc -l < ont.fasta.bam.bed.self.tsv` = 3117258  
TP = `comm -12 ont.kd.out.tsv ont.fasta.bam.bed.self.tsv | wc -l` = 1360726 (2615327)  
FP = `comm -23 ont.kd.out.tsv ont.fasta.bam.bed.self.tsv | wc -l` = 204942 (749711)  
FN = `comm -13 ont.kd.out.tsv ont.fasta.bam.bed.self.tsv | wc -l` = 1756532 (501931)  
Sensitivity: TP/P = TP/(TP+FN) = 43.7% (83.9%)  
Precision: TP/(TP+FP) = 86.9% (77.7%)  
F1 score: 2TP/(2TP+FP+FN) = 58.1% (80.7%)  

### PB
P = `wc -l < pacbio.fasta.bam.bed.self.tsv` = 11085073  
TP = `comm -12 pacbio.kd.out.tsv pacbio.fasta.bam.bed.self.tsv | wc -l` = 2179694 (5122592)  
FP = `comm -23 pacbio.kd.out.tsv pacbio.fasta.bam.bed.self.tsv | wc -l` = 540605 (788333)  
FN = `comm -13 pacbio.kd.out.tsv pacbio.fasta.bam.bed.self.tsv | wc -l` = 8905379 (5962481)  
Sensitivity: TP/P = TP/(TP+FN) = 19.7% (46.2%)  
Precision: TP/(TP+FP) = 80.1% (86.7%)  
F1 score: 2TP/(2TP+FP+FN) = 31.6% (60.3%)  

### Resulting table
Expanded Table 2 of Chu *et al.*, *Bioinformatics* 2017 (https://doi.org/10.1093/bioinformatics/btw811), rows *MHAP (default)* and *kd* added. Summary: kd is less sensitive (by design) but as precise.

|  | PB ||| ONT |||
|---|---|---|---|---|---|---|
|  | Sens. (%) | Prec. (%) | F1 (%) | Sens. (%) | Prec. (%) | F1 (%) |
| BLASR | 66.0 | 96.5 | 78.3 | 89.9 | 73.0 | 80.6 |
| DALIGNER | 83.8 | 85.8 | 84.8 | 92.9 | 91.0 | 91.9 |
| MHAP | 79.8 | 79.8 | 79.8 | 91.2 | 82.0 | 86.3 |
| **MHAP (default)** | **46.2** | **86.7** | **60.3** | **83.9** | **77.7** | **80.7** |
| GraphMap | 71.7 | 94.0 | 81.4 | 90.6 | 93.4 | 92.0 |
| Minimap | 59.6 | 83.8 | 69.7 | 91.2 | 95.4 | 93.2 |
| **kd** | **19.7** | **80.1** | **31.6** | **43.7** | **86.9** | **58.1** |
