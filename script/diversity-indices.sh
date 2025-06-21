# Generate diversity indices using the q2-diversity plugin in the QIIIME 2 implementation

# Generate alpha diversity (e.g., Shannon)
qiime diversity alpha \
  --i-table <table.qza> \
  --p-metric <shannon> \ #<ace> for the ACE index
  --o-alpha-diversity <shannon-diversity.qza>

# Generate Faith’s Phylogenetic Diversity (PD)
qiime diversity alpha-phylogenetic \
  --i-phylogeny <rooted-tree.qza> \
  --i-table <table.qza> \
  --p-metric <faith_pd> \
  --o-alpha-diversity <faith-pd-diversity.qza>

# Generate beta diversity distance matrix (e.g., Bray-Curtis)
qiime diversity beta \
  --i-table <feature-table.qza> \
  --p-metric braycurtis \
  --o-distance-matrix <braycurtis-distance.qza>

# Run PCoA on Bray-Curtis distance matrix
qiime diversity pcoa \
  --i-distance-matrix <braycurtis-distance.qza> \
  --o-pcoa <braycurtis-pcoa-results.qza>

# Visualize diversity group significance using ggplot2 in R
