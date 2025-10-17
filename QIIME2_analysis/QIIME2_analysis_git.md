# QIIME2 analysis

## Import data into QIIME2


```bash
## importing manually to slightly edit sample_ids due to repeats
```


```bash
#!/bin/bash

# load
module load miniconda3
conda activate qiime2-2022.8

#Script to generate manifest file from gzipped raw fastq files 

DATA_DIR="/srv/projects/microbial_inducers/data/all_reads_16S"
OUT_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S/01_import"

echo "Creating manifest file"

echo -e "sample-id,absolute-filepath,direction" > $OUT_DIR/manifest.txt #initiates manifest file and prints header

for i in $DATA_DIR/*R1_001.fastq.gz; do
    echo $i
    name=$(basename $i | cut -d '_' -f 1,3)
    echo $name
    echo "${name},${i},forward" >>  $OUT_DIR/manifest.txt #prints info and adds it to manifest
    echo "${name},${i/_R1/_R2},reverse" >>  $OUT_DIR/manifest.txt #prints print again but substitutes R1 for R2
done

echo "Done"
echo "Importing data into QIIME2"

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $OUT_DIR/manifest.txt \
--output-path $OUT_DIR/paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33

echo "Done"

conda deactivate
```

## Process 16S data using DADA2 pipeline

### Summarise sequence data


```bash
#!/bin/bash

# load
module load miniconda3
conda activate qiime2-2022.8

QIIME_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S/"

# summarise to get sequuence quality
qiime demux summarize \
  --i-data $QIIME_DIR/01_import/paired-end-demux.qza \
  --o-visualization $QIIME_DIR/02_get_ASVs/demux.qzv
conda deactivate
```

### Remove primer sequences


```bash
#!/bin/bash

# adapters and barcodes already removed
# removing only primer sequences
# load
module load miniconda3
conda activate qiime2-2022.8

QIIME_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S/"

# trim primer seqs
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences $QIIME_DIR/01_import/paired-end-demux.qza \
 --p-cores 16 \
 --p-front-f GTGYCAGCMGCCGCGGTAA \
 --p-front-r GGACTACNVGGGTWTCTAAT \
 --o-trimmed-sequences $QIIME_DIR/02_get_ASVs/paired-end-demux-trimmed.qza
 
# summarise to check if reads have been trimmed
qiime demux summarize \
  --i-data $QIIME_DIR/02_get_ASVs/paired-end-demux-trimmed.qza \
  --o-visualization $QIIME_DIR/02_get_ASVs/demux-trimmed.qzv
  
# add check read files
qiime tools export \
  --input-path $QIIME_DIR/02_get_ASVs/paired-end-demux-trimmed.qza \
  --output-path $QIIME_DIR/02_get_ASVs/paired-end-demux-trimmed
conda deactivate

```

### Run DADA2 to get ASVs


```bash
#!/bin/bash

# load
module load miniconda3
conda activate qiime2-2022.8

# running on mrca

QIIME_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S/02_get_ASVs"

# pick ASVs
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $QIIME_DIR/paired-end-demux-trimmed.qza \
  --p-trunc-len-f 235 \
  --p-trunc-len-r 203 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-n-threads 16 \
  --o-table $QIIME_DIR/table-dada2-3.qza \
  --o-representative-sequences $QIIME_DIR/rep-seqs-dada2-3.qza \
  --o-denoising-stats $QIIME_DIR/stats-dada2-3.qza

# summarise
qiime metadata tabulate \
  --m-input-file $QIIME_DIR/stats-dada2-3.qza \
  --o-visualization $QIIME_DIR/stats-dada2-3.qzv
conda deactivate

```
#mean
percentage of input passed filter =	89.41886076
percentage of input merged =	80.58850211
percentage of input non-chimeric =	80.16204641

### Summarise feature table 


```bash
#!/bin/bash

# load
module load miniconda3
conda activate qiime2-2022.8

# running on mrca

QIIME_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S/02_get_ASVs"

qiime feature-table summarize \
  --i-table $QIIME_DIR/table-dada2.qza \
  --o-visualization $QIIME_DIR/table-dada2.qzv \
  --m-sample-metadata-file $QIIME_DIR/metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data $QIIME_DIR/rep-seqs-dada2.qza \
  --o-visualization $QIIME_DIR/rep-seqs-dada2.qzv
```

# Get phylogeny

### Generate a tree for phylogenetic diversity analyses


```bash
#!/bin/bash

# load
module load miniconda3
conda activate qiime2-2022.8

# running on mrca

QIIME_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S"

# get phylogeny
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $QIIME_DIR/02_get_ASVs/rep-seqs-dada2.qza \
  --o-alignment $QIIME_DIR/03_phylogeny_and_diversity/aligned-rep-seqs.qza \
  --o-masked-alignment $QIIME_DIR/03_phylogeny_and_diversity/masked-aligned-rep-seqs.qza \
  --o-tree $QIIME_DIR/03_phylogeny_and_diversity/unrooted-tree.qza \
  --o-rooted-tree $QIIME_DIR/03_phylogeny_and_diversity/rooted-tree.qza
```

# Taxonomic analysis

### Assign taxonomy to 16S reads using the q2-feature-classifier


```bash
# reference database downloaded from https://docs.qiime2.org/2022.8/data-resources/
# include full length trained 16S Silva 138 SSU db
```

### Train and classify reads silva 138.1 release


```bash
#!/bin/bash

# load
module load miniconda3
conda activate qiime2-2022.8

# running on mrca

QIIME_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S/04_taxonomic_analysis"
SEQ_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S/02_get_ASVs"

#note: --p-min-length, --p-max-length are based on rep-seqs.qzv file
# qiime recommends not to truncate for paired reads due to the variable lengths

# extract reads from DB
qiime feature-classifier extract-reads \
  --i-sequences $QIIME_DIR/silva-138-99-full-trained-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
#  --p-trunc-len \
  --p-min-length 235 \
  --p-max-length 426 \
  --o-reads $QIIME_DIR/ref-seqs.qza

# train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $QIIME_DIR/ref-seqs.qza \
  --i-reference-taxonomy $QIIME_DIR/silva-138-99-full-trained-tax.qza \
  --o-classifier $QIIME_DIR/classifier.qza

# run classifier
qiime feature-classifier classify-sklearn \
  --i-classifier $QIIME_DIR/classifier.qza \
  --i-reads $SEQ_DIR/rep-seqs-dada2.qza \
  --o-classification $QIIME_DIR/taxonomy.qza

# summarise results
qiime metadata tabulate \
  --m-input-file $QIIME_DIR/taxonomy.qza \
  --o-visualization $QIIME_DIR/taxonomy.qzv
  
# plot taxa as bar plot
qiime taxa barplot \
  --i-table $SEQ_DIR/table-dada2.qza \
  --i-taxonomy $QIIME_DIR/taxonomy.qza \
  --m-metadata-file $QIIME_DIR/../01_import/metadata.tsv \
  --o-visualization $QIIME_DIR/taxa-bar-plots.qzv

conda deactivate
```

#### Plot taxa


```bash
#!/bin/bash

# load
module load miniconda3
conda activate qiime2-2022.8

# running on mrca

QIIME_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S/04_taxonomic_analysis"

# plot taxa as bar plot
qiime taxa barplot \
  --i-table $QIIME_DIR/02_get_ASV/table-dada2.qza \
  --i-taxonomy $QIIME_DIR/taxonomy.qza/taxonomy.qza \
  --m-metadata-file $QIIME_DIR/01_impor/metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

conda deactivate

```

# Export data for downstream analysis


```bash
#!/bin/bash

# load
module load miniconda3
conda activate qiime2-2022.8

QIIME_DIR="/srv/projects/microbial_inducers/analysis/qiime2_16S/"

qiime tools export \
    --input-path $QIIME_DIR/02_get_ASVs/table-dada2.qza  \
    --output-path $QIIME_DIR/05_exported_files/table-dada2
qiime tools export \
    --input-path $QIIME_DIR/02_get_ASVs/rep-seqs-dada2.qza \
    --output-path $QIIME_DIR/05_exported_files/rep-seqs
qiime tools export \
    --input-path $QIIME_DIR/03_phylogeny_and_diversity/rooted-tree.qza  \
    --output-path $QIIME_DIR/05_exported_files/rooted-tree
qiime tools export \
    --input-path $QIIME_DIR/03_phylogeny_and_diversity/unrooted-tree.qza  \
    --output-path $QIIME_DIR/05_exported_files/unrooted-tree
qiime tools export \
    --input-path $QIIME_DIR/04_taxonomic_analysis/taxonomy.qza  \
    --output-path $QIIME_DIR/05_exported_files/taxonomy

conda deactivate

# to convert biom tsv file (feature table)
biom convert \
-i $QIIME_DIR/05_exported_files/table-dada2/feature-table.biom \
-o $QIIME_DIR/05_exported_files/table-dada2/feature-table.tsv \
--to-tsv

```
