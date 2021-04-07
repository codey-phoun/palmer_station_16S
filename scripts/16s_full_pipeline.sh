#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=4500
#SBATCH --output=16s_qiime2_full_pipeline.log
#SBATCH --job-name=16s_full

module load intel-python3
source activate qiime2-2020.11

script_start=`date +%s`
job_start=`date +%s`

directory=/home/cphoun/palmer_station

# Remove adapters and primers from 16s fastq files with Cutadapt
# List of 16S fastq libraries stored in 16S_Libraries.txt
for library in `cat ${directory}/metadata/16S_Libraries.txt`
do
    echo "Processing $library"

    cutadapt \
    -b GTGCCAGCMGCCGCGGTAA `# FWD Read 1 515F (Caporaso)` \
    -b TTACCGCGGCKGCTGGCAC `# FRC` \
    -b GGACTACHVGGGTWTCTAAT `# REV 806R (Caporaso)` \
    -b ATTAGAWACCCBDGTAGTCC `# RRC` \
    -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA `# LTHT_1` \
    -b TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT `# LTHT_1_RC` \
    -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT `# LTHT_2` \
    -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT `# LTHT_2_RC` \
    -b TGGAATTCTCGGGTGCCAAGG `# SmallRNA` \
    -b CCTTGGCACCCGAGAATTCCA `# SmallRNA_RC` \
    -b CTGTCTCTTATACACATCT `# Nextera_1` \
    -b AGATGTGTATAAGAGACAG `# Nextera_2`\
    -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT `# TruSeq Universal Adapter` \
    -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG `# TruSeq Index Adapter` \
    -B GTGCCAGCMGCCGCGGTAA `# FWD Read 2` \
    -B TTACCGCGGCKGCTGGCAC `# FRC` \
    -B GGACTACHVGGGTWTCTAAT `# REV` \
    -B ATTAGAWACCCBDGTAGTCC `# RRC` \
    -B AGATCGGAAGAGCACACGTCTGAACTCCAGTCA `# LTHT_1` \
    -B TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT `# LTHT_1_RC` \
    -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT `# LTHT_2` \
    -B ACACTCTTTCCCTACACGACGCTCTTCCGATCT `# LTHT_2_RC` \
    -B TGGAATTCTCGGGTGCCAAGG `# SmallRNA` \
    -B CCTTGGCACCCGAGAATTCCA `# SmallRNA_RC` \
    -B CTGTCTCTTATACACATCT `# Nextera_1` \
    -B AGATGTGTATAAGAGACAG `# Nextera_2`\
    -B AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT `# TruSeq Universal Adapter` \
    -B GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG `# TruSeq Index Adapter` \
    --cores 0 `# auto-detect number of CPU cores to use` \
    --overlap 10 `# minimum length overlap between read and adapter` \
    --times 4 `# max number of adapters removed from each read` \
    --quality-cutoff 20,20 `# trim low-quality bases from 5' and 3' ends before adapter removal` \
    -o ${directory}/trimmed/${library}_1_trimmed.fastq `# Read 1 output` \
    -p ${directory}/trimmed/${library}_2_trimmed.fastq `# Read 2 output` \
    ${directory}/data/${library}_1.fastq `# Read 1 input` \
    ${directory}/data/${library}_2.fastq `# Read 2 input` \
    > ${directory}/trimmed/${library}_trim_log.txt
done

echo "Finished adapter and primer trimming with Cutadapt"
job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Cutadapt runtime: $hours:$minutes:$seconds (hh:mm:ss)"

# Import 16s libraries individually to QIIME2 for processing
# Trimmed 16s fastq filepaths need to be stored in the data_manifest_${library}.txt metadata file
echo "Importing 16s libraries to QIIME2"

for season in S1 S2 S3
do
    job_start=`date +%s`

    qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path ${directory}/metadata/data_manifest_${season}.txt \
    --output-path ${directory}/qiime2/${season}_PE_trimmed.qza \
    --input-format PairedEndFastqManifestPhred33V2

    qiime demux summarize \
    --i-data ${directory}/qiime2/${season}_PE_trimmed.qza \
    --p-n 1000000 \
    --o-visualization ${directory}/visualizations/${season}_summary.qzv

    echo "Finished importing ${season} into QIIME2"
    job_end=`date +%s`
    runtime=$((job_end-job_start))
    hours=$((runtime / 3600))
    minutes=$(( (runtime % 3600) / 60 ))
    seconds=$(( (runtime % 3600) % 60 ))
    echo "${season} import runtime: $hours:$minutes:$seconds (hh:mm:ss)"
done

# Process 16s libraries with DADA2 in QIIME2 to create ASVs
# Denoises paired-end sequences, dereplicates, and filters chimeras
echo "Processing 16s libraries with DADA2"

for season in S1 S2 S3
do
    job_start=`date +%s`

    qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ${directory}/qiime2/${season}_PE_trimmed.qza \
    --p-trim-left-f 0 `# do not truncate or trim reads, already processed with Cutadapt` \
    --p-trim-left-r 0 \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-max-ee-f 4 `# reads with expected errors higher than this are discarded` \
    --p-max-ee-r 4 \
    --p-min-fold-parent-over-abundance 8 `# minimum abundance of potential parents of chimera, increase value to reduce chimera detection sensitivity` \
    --p-n-threads 0 `# all available cores will be used` \
    --p-pooling-method 'pseudo' `# approximate pooling of sample`\
    --p-chimera-method 'consensus' `# default chimera detection method`\
    --p-n-reads-learn 2500000 `# number of reads to use when training the error model`\
    --o-table ${directory}/qiime2/${season}_feature_table.qza \
    --o-representative-sequences ${directory}/qiime2/${season}_rep_seqs.qza \
    --o-denoising-stats ${directory}/qiime2/${season}_denoising_stats.qza \
    --verbose \
    --p-no-hashed-feature-ids

    # view summary statistics of DADA2 run
    qiime metadata tabulate \
    --m-input-file ${directory}/qiime2/${season}_denoising_stats.qza \
    --o-visualization ${directory}/visualizations/${season}_denoising_stats.qzv \

    echo "Finished processing ${season} with DADA2"
    job_end=`date +%s`
    runtime=$((job_end-job_start))
    hours=$((runtime / 3600))
    minutes=$(( (runtime % 3600) / 60 ))
    seconds=$(( (runtime % 3600) % 60 ))
    echo "${season} DADA2 runtime: $hours:$minutes:$seconds (hh:mm:ss)"
done

# Merge all sequences and feature tables
echo "Merging 16s sequences and feature tables"
job_start=`date +%s`

library=16S_libraries

# merge all DADA2 sequences
qiime feature-table merge-seqs \
--i-data ${directory}/qiime2/S1_rep_seqs.qza \
--i-data ${directory}/qiime2/S2_rep_seqs.qza \
--i-data ${directory}/qiime2/S3_rep_seqs.qza \
--o-merged-data ${directory}/qiime2/${library}_rep_seqs.qza \

# merge all DADA2 feature tables
qiime feature-table merge \
--i-tables ${directory}/qiime2/S1_feature_table.qza \
--i-tables ${directory}/qiime2/S2_feature_table.qza \
--i-tables ${directory}/qiime2/S3_feature_table.qza \
--o-merged-table ${directory}/qiime2/${library}_feature_table.qza \

# information on how many sequences are associated with each sample and with each feature, 
# histograms of those distributions, and some related summary statistics
qiime feature-table summarize \
--i-table ${directory}/qiime2/${library}_feature_table.qza \
--o-visualization ${directory}/visualizations/${library}_feature_table.qzv \

# mapping of feature IDs to sequences, and provide links to 
# easily BLAST each sequence against the NCBI nt database
qiime feature-table tabulate-seqs \
--i-data ${directory}/qiime2/${library}_rep_seqs.qza \
--o-visualization ${directory}/visualizations/${library}_rep_seqs.qzv \

job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Merge runtime: $hours:$minutes:$seconds (hh:mm:ss)"

# Classify ASVs with VSEARCH and silva-138
echo "Performing taxanomic classification with VSEARCH and silva-138"
job_start=`date +%s`

qiime feature-classifier classify-consensus-vsearch \
--i-query ${directory}/qiime2/${library}_rep_seqs.qza \
--i-reference-reads ${directory}/metadata/silva-138-ssu-nr99-dna-seqs.qza \
--i-reference-taxonomy ${directory}/metadata/silva-138-ssu-nr99-tax.qza \
--p-maxaccepts 1000 \
--p-top-hits-only \
--p-perc-identity 0.8 \
--p-query-cov 0.8 \
--p-strand "both" \
--p-threads 56 \
--o-classification ${directory}/qiime2/${library}_vsearch_taxonomy.qza

# Visualize taxonomy results table
qiime metadata tabulate \
--m-input-file ${directory}/qiime2/${library}_vsearch_taxonomy.qza \
--o-visualization ${directory}/visualizations/${library}_vsearch_taxonomy.qzv

# Merge feature_table and taxonomy results
qiime feature-table transpose \
--i-table ${directory}/qiime2/${library}_feature_table.qza \
--o-transposed-feature-table ${directory}/qiime2/${library}_feature_table_transposed.qza 

qiime metadata tabulate \
--m-input-file ${directory}/qiime2/${library}_feature_table_transposed.qza \
--m-input-file ${directory}/qiime2/${library}_vsearch_taxonomy.qza \
--o-visualization ${directory}/visualizations/${library}_combined_taxonomy_feature_table.qzv

# Filter out mitochondria, chloroplast, eukaryota, and unassigned results:
qiime taxa filter-seqs \
--i-sequences ${directory}/qiime2/${library}_rep_seqs.qza \
--i-taxonomy ${directory}/qiime2/${library}_vsearch_taxonomy.qza \
--p-include d__ \
--p-exclude mitochondria,chloroplast,eukaryota \
--o-filtered-sequences ${directory}/qiime2/${library}_sequences_clean.qza

qiime taxa filter-table \
--i-table ${directory}/qiime2/${library}_feature_table.qza \
--i-taxonomy ${directory}/qiime2/${library}_vsearch_taxonomy.qza \
--p-include d__ \
--p-exclude mitochondria,chloroplast,eukaryota \
--o-filtered-table ${directory}/qiime2/${library}_feature_table_clean.qza

## View filtered results
qiime feature-table tabulate-seqs \
--i-data ${directory}/qiime2/${library}_sequences_clean.qza \
--o-visualization ${directory}/visualizations/${library}_sequences_clean.qzv

qiime feature-table summarize \
--i-table ${directory}/qiime2/${library}_feature_table_clean.qza \
--o-visualization ${directory}/visualizations/${library}_feature_table_clean.qzv

# Merge filtered feature_table and taxonomy results
qiime feature-table transpose \
--i-table ${directory}/qiime2/${library}_feature_table_clean.qza \
--o-transposed-feature-table ${directory}/qiime2/${library}_feature_table_clean_transposed.qza 

qiime metadata tabulate \
--m-input-file ${directory}/qiime2/${library}_feature_table_clean_transposed.qza  \
--m-input-file ${directory}/qiime2/${library}_vsearch_taxonomy.qza \
--o-visualization ${directory}/visualizations/${library}_combined_taxonomy_feature_table_clean.qzv

## View taxonomic composition
qiime taxa barplot \
--i-table ${directory}/qiime2/${library}_feature_table_clean.qza \
--i-taxonomy ${directory}/qiime2/${library}_vsearch_taxonomy.qza \
--m-metadata-file ${directory}/metadata/16S_metadata.tsv \
--o-visualization ${directory}/visualizations/${library}_taxa_bar_plots.qzv

job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Taxonomic classification runtime: $hours:$minutes:$seconds (hh:mm:ss)"

# Create phylogenetic tree for diversity analysis
echo "Creating phylogenetic tree"
job_start=`date +%s`

qiime phylogeny align-to-tree-mafft-iqtree \
--i-sequences ${directory}/qiime2/${library}_sequences_clean.qza \
--p-n-threads 'auto' \
--o-alignment ${directory}/phylogeny/${library}_aligned_sequences.qza \
--o-masked-alignment ${directory}/phylogeny/${library}_masked_aligned_sequences.qza \
--o-tree ${directory}/phylogeny/${library}_iqtree.qza \
--o-rooted-tree ${directory}/phylogeny/${library}_iqtree_rooted.qza \

job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Tree construction with iqtree runtime: $hours:$minutes:$seconds (hh:mm:ss)"

# Perform core diversity metrics
echo "Calculating core diversity metrics"
job_start=`date +%s`

qiime diversity core-metrics-phylogenetic \
--i-phylogeny ${directory}/phylogeny/${library}_iqtree_rooted.qza \
--i-table ${directory}/qiime2/${library}_feature_table_clean.qza \
--p-sampling-depth 193308 \
--m-metadata-file ${directory}/metadata/16S_metadata.tsv \
--output-dir ${directory}/${library}_diversity_results/ \

job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Diversity calculations runtime: $hours:$minutes:$seconds (hh:mm:ss)"

# Perform rarefaction plotting
echo "Plotting rarefaction curves"
job_start=`date +%s`

qiime diversity alpha-rarefaction \
--i-table ${directory}/qiime2/${library}_feature_table_clean.qza \
--i-phylogeny ${directory}/phylogeny/${library}_iqtree_rooted.qza \
--p-max-depth 594964 `# largest library size from clean feature table`\
--m-metadata-file ${directory}/metadata/16S_metadata.tsv \
--o-visualization ${directory}/${library}_diversity_results/alpha-rarefaction.qzv \

job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Rarefaction calculations runtime: $hours:$minutes:$seconds (hh:mm:ss)"

script_end=`date +%s`
runtime=$((script_end-script_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "16s full pipeline runtime: $hours:$minutes:$seconds (hh:mm:ss)"