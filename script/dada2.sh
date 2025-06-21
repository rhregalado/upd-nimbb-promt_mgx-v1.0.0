# DADA2 script using QIIME 2 implementation
# Install QIIME 2 in your system. Documentation can be found here: https://docs.qiime2.org/2024.10/

# Activate QIIME 2

# Import paired-end fastq files
# Ensure the manifest file is in your directory
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path <manifest.csv> \
  --output-path <paired-end-demux.qza> \
  --input-format PairedEndFastqManifestPhred33
  
# Summarize imported data
qiime demux summarize \
  --i-data <paired-end-demux.qza> \
  --o-visualization <paired-end-demux.qzv>
  
# Quality filtering
qiime quality-filter q-score \
  --i-demux <paired-end-demux.qza> \
  --o-filtered-sequences <demux-filtered.qza> \
  --o-filter-stats <demux-filter-stats.qza>
  
# DADA2 denoising
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs <paired-end-demux.qza> \
  --p-trim-left-f <value> \
  --p-trim-left-r <value> \
  --p-trunc-len-f <value> \
  --p-trunc-len-r <value> \
  --o-table <table.qza> \
  --o-representative-sequences <rep-seqs.qza> \
  --o-denoising-stats <stats-dada2.qza>

# Tabulate denoising stats
qiime metadata tabulate \
  --m-input-file <stats-dada2.qza> \
  --o-visualization <stats-dada2.qzv>

# Feature table summary
qiime feature-table summarize \
  --i-table <table.qza> \
  --m-sample-metadata-file <metadata.tsv> \
  --o-visualization <table.qzv>

# Tabulate representative sequences
qiime feature-table tabulate-seqs \
  --i-data <rep-seqs.qza> \
  --o-visualization <rep-seqs.qzv>

# Generate phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences <rep-seqs.qza> \
  --o-alignment <aligned-rep-seqs.qza> \
  --o-masked-alignment <masked-aligned-rep-seqs.qza> \
  --o-tree <unrooted-tree.qza> \
  --o-rooted-tree <rooted-tree.qza>

# Alpha rarefaction
qiime diversity alpha-rarefaction \
  --i-table <table.qza> \
  --i-phylogeny <rooted-tree.qza> \
  --p-max-depth <value> \
  --m-metadata-file <metadata.tsv> \
  --o-visualization <alpha-rarefaction.qzv>

# Core diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny <rooted-tree.qza> \
  --i-table <table.qza> \
  --p-sampling-depth <value> \
  --m-metadata-file <metadata.tsv> \
  --output-dir <core-metrics-results>

# Taxonomic classification (16S). Extract classifier from the mapfiles folder
qiime feature-classifier classify-sklearn \
  --i-classifier silva-classifier.qza \
  --i-reads <rep-seqs.qza> \
  --o-classification <taxonomy.qza>

# Tabulate taxonomy
qiime metadata tabulate \
  --m-input-file <taxonomy.qza> \
  --o-visualization <taxonomy.qzv>

# Taxonomy barplot
qiime taxa barplot \
  --i-table <table.qza> \
  --i-taxonomy <taxonomy.qza> \
  --m-metadata-file <metadata.tsv> \
  --o-visualization <taxa-bar-plots.qzv>
