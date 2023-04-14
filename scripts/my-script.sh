#!/usr/bin/env bash
# make empty directories in git repository. Theses are the analysis, DNAseqAnalysis, docs and data directories.
mkdir -p DNAseqAnalysis analysis docs data
# add a README.md to each directory
# the scripts directory already exists
for my_directory in scripts DNAseqAnalysis analysis docs data;do
  touch ${my_directory}/README.md
  echo "# ${my_directory}" >> ${my_directory}/README.md
done
