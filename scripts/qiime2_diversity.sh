#! /bin/bash
#
#$ -S /bin/bash
#$ -r y
#$ -j y
#$ -R yes


input="Deblur_all_feature-table-1152samples_10reads_filtered_rep_seqs.biom"
group=$(sed s/.biom//g <<< ${input})
SCRATCH_JOB=/scratch/$USER/jobs/$JOB_ID
mkdir -p $SCRATCH_JOB


cp -r ${input} $SCRATCH_JOB/
cd $SCRATCH_JOB
mkdir ${group}

qiime tools import \
--type 'FeatureTable[Frequency]' \
--input-path  Deblur_all_feature-table-1152samples_10reads_filtered_rep_seqs.biom \
--input-format BIOMV210Format \
--output-path otu_table_filtered.qza

qiime tools import --input-path rep_seqs.fasta --output-path rep_sequences.qza --type 'FeatureData[Sequence]'

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep_sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table otu_table_filtered.qza \
  --p-sampling-depth 10000 \
  --m-metadata-file map.txt \
  --output-dir core-metrics-results
  
  # alpha rarefaction
  qiime diversity alpha-rarefaction \
  --i-table otu_table_filtered.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 10000 \
  --m-metadata-file map.txt \
  --o-visualization alpha-rarefaction.qzv
# rarefy feature table

qiime feature-table rarefy \
--i-table otu_table_filtered.qza \
--p-sampling-depth 10000 \
--o-rarefied-table otu_table_filtered10k.qza

# alpha diversity
qiime diversity alpha \
  --i-table otu_table_filtered.qza \
  --p-metric shannon \
  --o-alpha-diversity shannon_vector.qza
  
  qiime diversity alpha \
  --i-table otu_table_filtered.qza \
  --p-metric chao1 \
  --o-alpha-diversity chao1_vector.qza
  
  
# extract the files from the qzv
for FILE in $(ls *.qzv); do unzip $FILE -d ${FILE%%.*}; done
for FILE in $(ls *.qza); do unzip $FILE -d ${FILE%%.*}; done

mv * ${group}
mv -f ${group} /home/diversity/
 

