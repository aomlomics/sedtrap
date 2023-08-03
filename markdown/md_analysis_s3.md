### Goal: Run Tourmaline on GMT-1 16s and 18s data to generate figures and tables for further analysis in R. 

Platform: 
This code was run on the NOAA parallel works AWS cloud instance with the following parameters- 
* Controller node: m6id.2xlarge (8 vCPUs, 32 GB Memory, amd64) 
* Image disk settings: 1 image disk, 250 GB
* Compute node: c3.8xlarge (32 vCPUs, 60 GB Memory, amd64) max nodes, 2

### Notes:
* Reference database used: silva-138-99-seqs-515-806

## Clone directory
```
cd /contrib/Rachael.Storo/16S-GMT1/tourmaline-GMT-1_16S
git clone https://github.com/aomlomics/tourmaline.git

mv tourmaline tourmaline-GMT-1_16S
```

## Dependencies

```
wget https://data.qiime2.org/distro/core/qiime2-2023.2-py38-osx-conda.yml
conda env create -n qiime2-2023.2 --file qqiime2-2023.2-py38-osx-conda.yml
```
```
conda activate qiime2-2023.2
conda install -c conda-forge -c bioconda snakemake biopython muscle clustalo tabulate
conda install -c conda-forge deicode
pip install empress
qiime dev refresh-cache
conda install -c bioconda bioconductor-msa bioconductor-odseq
```

## Initiate Tourmaline
```
cd /contrib/Rachael.Storo/16S-GMT1/tourmaline-GMT-1_16S01-imported
```
```
wget https://data.qiime2.org/2021.2/common/silva-138-99-seqs-515-806.qza
wget https://data.qiime2.org/2021.2/common/silva-138-99-tax-515-806.qza
```
```
mv silva-138-99-tax-515-806.qza reftax.qza
mv silva-138-99-seqs-515-806.qza refseqs.qza
```
```
cd /contrib/Rachael.Storo/16S-GMT1/tourmaline-GMT-1_16S
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
cd /contrib/Rachael.Storo/16S-GMT1/tourmaline-GMT-1_16S
snakemake dada2_se_denoise --cores all
```
```
snakemake dada2_se_taxonomy_unfiltered --cores all
snakemake dada2_se_taxonomy_filtered --cores all
```
```
snakemake dada2_se_diversity_unfiltered --cores all
snakemake dada2_se_diversity_filtered --cores all
```
```
snakemake dada2_pe_report_unfiltered
```












