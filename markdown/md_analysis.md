### Goal: Run Tourmaline on GMT-1 16s and 18s data to generate figures and tables for further analysis in R. 

Platform: 
This code was run on an intel mac with the following parameters- 
* 2.6 GHz 6-core intel core i7
* AMD Radeon Pro 5300M, 4GB, Intel UHD Graphics 630 1536 MB
* 16 GB 2667 MHz DDR4
* Ventura 13.2.1

### Notes:
* 16S Reference database used: silva-183_1-99-515f_926r-sediment-saline-classifier.qza
* 18S Reference database used: (18s) PR2 v. 5.0.0

## Clone directory
```
git clone https://github.com/aomlomics/tourmaline.git
mv tourmaline tourmaline-GMT-1_16S 
```
## Dependencies
```
conda create -c conda-forge -c bioconda -n snakemake snakemake-minimal

wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-osx-conda.yml
conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-osx-conda.yml

conda activate qiime2-2023.5
conda install -c conda-forge -c bioconda biopython muscle clustalo tabulate
conda install -c conda-forge deicode
pip install empress
qiime dev refresh-cache
conda install -c bioconda bioconductor-msa bioconductor-odseq
```
## Reference Databases
### 18S
```
qiime tools import \
	--type 'FeatureData[Taxonomy]' \
	--input-path /Users/rachaelkarns/Downloads/pr2_version_5.0.0_SSU_mothur.tax \
	--output-path ref-taxonomy.qza \
	--input-format HeaderlessTSVTaxonomyFormat

qiime tools import \
	--type 'FeatureData[Sequence]' \
	--input-path /Users/rachaelkarns/Downloads/pr2_version_5.0.0_SSU_mothur.fasta \
	--output-path ref-sequence.qza

qiime feature-classifier \
	fit-classifier-naive-bayes \
	--i-reference-reads ref-sequence.qza \
	--i-reference-taxonomy ref-taxonomy.qza \
	--o-classifier pr2-classifier.qza

qiime tools export \
	--input-path pr2-classifier.qza \
	--output-path refseqs.qza

mv refseqs.qza ../01-imported
```
## Trim Primers
### Trim primers- 18S
``` 
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ./ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end-18S-GMT-1.qza


qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end-18S-GMT-1.qza --p-adapter-f AGTAGGTGAACCTGCAGAAGGATC --p-adapter-r GACGGGCGGTGTGTAC --p-match-read-wildcards --p-match-adapter-wildcards --verbose --o-trimmed-sequences trimmed_remove_primers_wild-18S-GMT-1.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences trimmed_remove_primers_wild-18S-GMT-1.qza --p-front-f GTACACACCGCCCGTC --p-front-r TGATCCTTCTGCAGGTTCACCTAC --p-match-read-wildcards --p-match-adapter-wildcards --p-discard-untrimmed --verbose --o-trimmed-sequences demux-trimmed-18S-GMT-1.qza

qiime tools export --input-path demux-trimmed-18S-GMT-1.qza --output-path trimmed-reads
```
### Trim primers- 16S
``` 
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ./ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end-16S-GMT-1.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end-16S-GMT-1.qza --p-anywhere-f GTGYCAGCMGCCGCGGTAA --p-anywhere-r CCGYCAATTYMTTTRAGTTT --p-match-read-wildcards TRUE --p-match-adapter-wildcards TRUE --verbose --o-trimmed-sequences trimmed_remove_primers_wild-16S-GMT-1.qza

qiime tools export --input-path trimmed_remove_primers_wild-16S-GMT-1.qza --output-path trimmed-reads
```
## Run Tourmaline

```
cd ../tourmaline-GMT-1_16S
snakemake --use-conda dada2_se_denoise --cores all
snakemake --use-conda dada2_se_taxonomy_unfiltered --cores all
snakemake --use-conda dada2_se_taxonomy_filtered --cores all
snakemake --use-conda dada2_se_diversity_unfiltered --cores all
snakemake --use-conda dada2_se_diversity_filtered --cores all
snakemake --use-conda dada2_pe_report_unfiltered

cd ../tourmaline-GMT-1_18S
snakemake --use-conda dada2_se_denoise --cores all
snakemake --use-conda dada2_se_taxonomy_unfiltered --cores all
snakemake --use-conda dada2_se_taxonomy_filtered --cores all
snakemake --use-conda dada2_se_diversity_unfiltered --cores all
snakemake --use-conda dada2_se_diversity_filtered --cores all
snakemake --use-conda dada2_pe_report_unfiltered
```











