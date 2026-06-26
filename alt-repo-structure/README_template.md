# Sediment Trap eDNA Time Series: Northern Gulf of Mexico

Analysis code and reproducible workflows for characterizing the biodiversity and temporal dynamics of the biological carbon pump using environmental DNA (eDNA) metabarcoding, applied to a five-year moored sediment trap time series in the northern Gulf of Mexico.

## Overview

This repository contains code for analyzing eDNA metabarcoding data (16S, 18S, 12S, COI rRNA genes) collected from dual moored sediment traps at two depths (≈550 m and ≈600 m) in the northern Gulf of Mexico (28°N, 89°W). Analyses integrate eDNA taxonomic profiles with complementary biogeochemical measurements (particulate flux, POC, PIC, opal), bio-optical backscattering data (b_bp), and discrete water column carbonate chemistry to characterize microbial and planktonic community composition driving carbon export.

**Core research questions:**
- How does the microbial and metazoan community composition of sinking particulate organic matter (POM) vary seasonally and interannually?
- Which taxa are consistently associated with carbon flux at mesopelagic depths, and can eDNA profiles serve as functional indicators of carbon export efficiency?
- How do eDNA-based community profiles compare with morphological (foraminiferal) and bio-optical (backscatter) approaches to characterizing the biological carbon pump (BCP)?
- Can we link eDNA community composition to ocean acidification (OA) and carbonate system variability?

## Dataset Description

**Collection & preservation:**
- Deployment: October 2019 – May 2025 (59 months, with annual mooring recovery/redeployment)
- Platform: Dual McLane PARFLUX Mark 78 automated sediment traps
- Depths: ≈120 m (Trap 1, euphotic zone; eDNA-optimized preservation) and ≈550–700 m (Trap 2, mesopelagic; formalin-fixed)
- Sampling interval: 2 weeks to 1 month (biweekly–monthly cup rotation)
- Total samples processed: 59 eDNA extractions from trap cup subsamples
- Coordinates: 28°N, 89–90°W (northern Gulf of Mexico, ~1,100–1,200 m water depth)

**Molecular data:**
- Markers: 16S rRNA (bacteria/archaea), 18S rRNA (eukaryotes), 12S rRNA (fish), COI (foraminifera and metazoans)
- Sequencing platform: Illumina (processed at Michigan State University RTSF)
- Processing: QIIME 2 / DADA2 denoising to amplicon sequence variants (ASVs)
- Taxonomy: Consensus-BLAST assignment against SILVA (16S/18S) and MetaZooGene (12S/COI) with percent-identity thresholds (≥95% species, 80--94% genus/LCA, <80% unassigned)
- Raw data: NCBI BioProject [INSERT], SRA Accessions [INSERT]

**Biogeochemical data:**
- Mass flux (mg·m²·d⁻¹)
- Particulate organic carbon (POC), nitrogen (PON), phosphorus (POP)
- Particulate inorganic carbon (PIC; carbonate)
- Biogenic silica (opal)
- Processing: Elemental analyzer (POC/PON), coulometry (PIC), spectrophotometry (opal)
- Provider: University of South Carolina (Claudia Benitez-Nelson lab)

**Bio-optical & carbonate chemistry:**
- Particulate backscattering (b_bp700) from moored sensor (QUA Real-Time Systems)
- CTD (T, S, depth), dissolved O₂ (Winkler)
- Discrete DIC, total alkalinity (TA), pH measurements
- Data archived: NCEI OCADS (oceanacidification.noaa.gov), NASA SeaBASS

**Metadata:**
- All samples follow FAIRe eDNA metadata checklist (v1.0.2)
- NCBI BioSample standardized attributes (host, environment, lat/lon, collection date)
- Sample-level QA/QC flags (extraction yield, PCR success, sequencing depth)

## Data Access & Archiving

**Public data repositories:**
- **Raw sequences:** NCBI Sequence Read Archive (SRA, BioProject [ID])
- **Taxonomy data:** OBIS/GBIF via edna2obis v2.0 pipeline (Darwin Core format)
- **Biogeochemistry & metadata:** NOAA National Centers for Environmental Information (NCEI) OADS
- **Code & workflows:** GitHub (this repository)

**Accessing processed data locally:**
```bash
# Clone repository
git clone https://github.com/[your-username]/sedtrap-edna-analysis.git
cd sedtrap-edna-analysis

# Download processed QIIME 2 artifacts and metadata
cd data/processed
# Instructions for downloading from Zenodo or OSF (if pre-archived)
bash download_data.sh
```

## Repository Structure

See `data/raw/README.md` for detailed data inventory and provenance.

```
scripts/01_qc_processing/       → Quality control, DADA2 ASV generation, taxonomy assignment
scripts/02_exploratory/         → Alpha/beta diversity, temporal trends, community composition
scripts/03_statistical_models/  → Seasonal GAMs, biogeochem correlation, eRNA analysis
scripts/04_multivariate/        → Integration with bio-optical data, BCP models
scripts/05_figures_tables/      → Manuscript figures and summary tables

data/processed/                 → QIIME 2 artifacts, CSV tables, metadata
results/figures/                → Generated manuscript & supplementary figures
results/tables/                 → Summary statistics and model outputs
results/statistics/             → Statistical test results, permanova outputs

notebooks/                      → R Markdown exploratory analyses & narrative
manuscripts/                    → Manuscript drafts and review response documents
docs/                           → Methods, SOP documentation, reproducibility guide
config/                         → Tourmaline, DADA2, and analysis configuration files
```

## Installation & Dependencies

**Requirements:**
- R ≥ 4.0 (phyloseq, vegan, tidyverse, ggplot2, ggvegan, gridExtra, mgcv)
- Python 3.8+ with QIIME 2 2023.7+ (via conda)
- Optional: Tourmaline workflow (Snakemake, Docker, or conda)

**Quick setup:**
```bash
# Create conda environment from YAML
conda env create -f environment.yml
conda activate sedtrap-edna

# Install R dependencies (from R console)
install.packages(c("phyloseq", "vegan", "tidyverse", "ggplot2", "gridExtra", "mgcv"))
```

**For full Tourmaline workflow:**
```bash
# See https://github.com/aomlomics/tourmaline for installation
# Or run via Docker: docker pull aomlomics/tourmaline:latest
```

## Quick Start

### 1. QC and ASV Generation (if running from raw FASTQ)
```bash
# Skip if using pre-processed QIIME 2 artifacts
cd scripts/01_qc_and_processing/

# Optimize DADA2 trimming parameters
Rscript trim_params_optimization.R \
  --metadata data/processed/metadata/sample_metadata.csv \
  --input_dir /path/to/raw_fastq

# Run full Tourmaline workflow
bash run_tourmaline.sh config/tourmaline_config.yaml
```

### 2. Exploratory Analysis
```bash
cd scripts/02_exploratory_analysis/

# Generate alpha and beta diversity metrics
Rscript alpha_diversity.R \
  --asv_table data/processed/asv_table.csv \
  --metadata data/processed/metadata/sample_metadata.csv \
  --output results/diversity_metrics.csv

Rscript beta_diversity_ordination.R \
  --asv_table data/processed/asv_table.csv \
  --metadata data/processed/metadata/sample_metadata.csv \
  --output results/figures/beta_diversity_nmds.pdf
```

### 3. Temporal Trends & Community Dynamics
```bash
# Analyze seasonal patterns and interannual variability
Rscript temporal_trends.R \
  --asv_table data/processed/asv_table.csv \
  --metadata data/processed/metadata/sample_metadata.csv \
  --tax_table data/processed/taxonomy_table.csv \
  --output results/temporal_summary.txt
```

### 4. Integration with Biogeochemistry
```bash
# Correlate eDNA profiles with flux and POC measurements
Rscript eDNA_vs_biogeochem.R \
  --asv_table data/processed/asv_table.csv \
  --flux_data data/processed/metadata/biogeochem_flux.csv \
  --output results/statistics/correlation_matrix.csv
```

### 5. Generate Manuscript Figures
```bash
cd scripts/05_figures_tables/

# Reproduce all main and supplementary figures
Rscript manuscript_fig1.R     # Time series + diversity
Rscript manuscript_fig2.R     # Composition heatmaps
Rscript manuscript_fig3.R     # Ordinations + drivers
Rscript supplementary_figures.R
```

## Key Analyses

### eDNA Metabarcoding Workflow
1. **Quality control:** Read filtering, chimera detection (DADA2)
2. **Denoising:** ASV inference (DADA2 algorithm; no OTU clustering)
3. **Taxonomy assignment:** Consensus-BLAST against SILVA (16S/18S) or MetaZooGene (12S/COI)
   - Species-level: ≥95% sequence identity
   - Genus/LCA: 80–94%
   - Unassigned: <80%
4. **Filtering:** Removal of rare ASVs (prevalence, min abundance)
5. **Diversity metrics:** Shannon index, Chao1 (alpha); Bray-Curtis (beta)

### Statistical Approaches
- **Seasonal decomposition:** Generalized additive models (GAM) with cyclic cubic splines
- **Multivariate analysis:** NMDS, Procrustes rotation (trap depths vs. time)
- **Biogeochemical correlation:** Spearman rank, Pearson (with significance testing)
- **Case studies:** 2022 Sargassum bloom (surface-to-depth export), foraminiferal eDNA validation

### Integration with Complementary Data
- **Bio-optical:** Cross-calibration of sediment trap b_bp700 with BGC-Argo float array
- **Carbonate chemistry:** Linking OA indices (Ω_arag, pH) to community composition
- **Flow cytometry:** Comparison with autonomous cell counts (where available)

## Publications

**Manuscripts in preparation:**
1. Anderson *et al.* – 3.5-year eDNA metabarcoding time series characterizing BCP community dynamics and seasonal export pulses (target: *Applied and Environmental Microbiology*)
2. Osborne *et al.* – Record 2022 Sargassum bloom: mesopelagic export and organic carbon pathways (target: *Nature Communications*)

**Related publications:**
- Thompson, L.R., *et al.* (2023). Decoding dissolved information: environmental DNA sequencing at global scale to monitor a changing ocean. *Curr. Opin. Biotechnol.* 81:102936.
- Silliman, K., Anderson, S., *et al.* (2023). A case study in sharing marine eDNA metabarcoding data to OBIS. *Biodivers. Inf. Sci. Stand.* 7:e111048.

## Contributing

Contributions, bug reports, and suggestions are welcome. Please:
1. Open an issue to describe the problem or feature
2. Create a branch for your work: `git checkout -b feature/your-idea`
3. Submit a pull request with a clear description
4. Follow existing code style and documentation conventions

For major changes, please open an issue first to discuss.

## Team

**Primary Investigators:**
- Luke R. Thompson (Mississippi State University / NOAA AOML)
- Emily Osborne (South Carolina Sea Grant / formerly NOAA AOML)
- Enrique Montes (CIMAS / NOAA AOML)

**Contributors:**
- Katherine Silliman (Research Scientist, eDNA protocols & data standards)
- Sean Anderson (formerly postdoc; eDNA extraction & bioinformatics)
- Rebecca Trinh (Research Scientist II; current analysis & manuscript preparation)
- Julie Richey (USGS; biogeochemistry & flux data)
- Claudia Benitez-Nelson (University of South Carolina; bulk sediment analysis)

**Funding:** NOAA Ocean Exploration and Research (OER), NOAA Ocean Acidification Program (OAP), NOAA Global Ocean Monitoring and Observing (GOMO), NOAA RESTORE Science Program, Mississippi State University

## License

This code is released under the [Apache 2.0 License](LICENSE). Please see LICENSE file for details.

Data are released under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/), with the exception of any data embargoed by NOAA or partner institutions (see data repository documentation).

## Citation

If using this repository or associated analyses, please cite:

Thompson, L.R., Anderson, S.R., Montes, E., *et al.* (2026). [Manuscript title]. *[Journal]*, [Volume/Issue]. doi: [INSERT]

For code, please cite the repository:

```bibtex
@software{thompson2026sedtrap,
  author = {Thompson, Luke R.},
  title = {sedtrap-edna-analysis: Sediment Trap eDNA Time Series Analysis},
  year = {2026},
  url = {https://github.com/[username]/sedtrap-edna-analysis},
  doi = {[INSERT]}
}
```

## Questions & Support

For questions about the repository:
- Open a GitHub issue for bug reports or feature requests
- Contact: luke.thompson@noaa.gov

For questions about the raw data or sample metadata:
- See `data/raw/README.md` for data access instructions
- NCBI BioProject: [INSERT]; SRA Accessions: [INSERT]

## Acknowledgments

We thank the captains, crews, and technicians of R/V Pelican, NOAA/OME, and the University of South Carolina for sediment trap deployment, recovery, and biogeochemical processing. Sequencing was performed at Michigan State University RTSF. This research was supported by NOAA OER (grant NA21OAR4320190-T3-01S006), NOAA OAP (grant [INSERT]), and the Mississippi State University Northern Gulf Institute. Emily Osborne was supported by South Carolina Sea Grant.
