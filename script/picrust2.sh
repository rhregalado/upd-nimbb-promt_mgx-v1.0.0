# Install PICRUSt2 in your system. Documentation can be found here: https://github.com/picrust/picrust2/wiki

# Activate PICRUSt2 in terminal

# Run PICRUSt2 pipeline
picrust2_pipeline.py \
  -s <ASVs-file.fna> \
  -i <ASVs-counts.biom> \
  -o <output-directory> \
  -p <number-of-threads>

# Add descriptions for each functional category 
# For Metacyc pathway abundance tables
add_descriptions.py \
  -i <path-to-pathways-abun-file.tsv.gz> \
  -m METACYC \
  -o <output-directory>/path_abun_unstrat_descrip.tsv.gz

# For KO pathway abundance tables
add_descriptions.py \
  -i <path-to-KO-metagenome-file.tsv.gz> \
  -m KO \
  -o <output-directory>/pred_metagenome_unstrat_descrip.tsv.gz

# For E.C. number pathway abundance tables (optional)
add_descriptions.py \
  -i <path-to-EC-metagenome-file.tsv.gz> \
  -m EC \
  -o <output-directory>/pred_metagenome_unstrat_descrip.tsv.gz
