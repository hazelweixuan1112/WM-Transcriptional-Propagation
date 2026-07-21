# WM-Transcriptional-Propagation

This repository contains analysis code for the manuscript **Correlated Gene Expression Supports Signal Propagation in Human Working Memory Networks**.

The analyses link working-memory (WM) activation, regional gene expression, intracranial stimulation-evoked response probability, receptor and macroscale brain maps, cell-type annotations, and developmental expression data to examine the molecular organization of signal propagation among WM-related brain regions.

## Repository Contents

The uploaded code is organized as a sequential analysis workflow:

| File | Purpose |
| --- | --- |
| `01.gene_selection_activate.ipynb` | Preprocesses AHBA expression and WM activation data, computes WM-related gene scores, selects top-ranked WM-matched genes, and prepares regional expression maps. |
| `02.gene_coexprssion.ipynb` | Constructs interregional gene co-expression / transcriptional similarity matrices from selected WM-related genes. |
| `03.corr.ipynb` | Aligns transcriptional similarity matrices with F-TRACT response probability matrices and computes their association. |
| `03-1.control.ipynb` | Performs control analyses, including spatial autocorrelation-preserving BrainSMASH null models and random gene controls. |
| `04.virtual_knockout.ipynb` | Performs virtual gene knockout analysis and calculates the Gene Contribution Index (GCI) for each selected gene. |
| `05.new_gene_score.ipynb` | Matches GCI gene groups with gene scores and prepares GCI-specific gene-score and regional-expression outputs. |
| `06.cell_type_deconvolution.ipynb` | Tests cell-type enrichment for GCI+ and GCI- gene groups using marker gene lists and permutation-based null models. |
| `07.PLS_posi_nega.ipynb` | Relates GCI gene-expression profiles to receptor, transporter, and macroscale brain maps using partial least squares (PLS) analysis. |
| `08.brainspan.ipynb` | Examines developmental expression patterns of GCI gene groups using BrainSpan data and generates developmental expression plots. |
| `09.PCA.ipynb` | Builds region-by-developmental-period expression matrices and performs PCA on developmental expression profiles. |
| `10.enrichment.R` | Performs GO, KEGG, and Reactome enrichment analyses for GCI+ and GCI- gene sets using `clusterProfiler` and `ReactomePA`. |

Files beginning with `._` are macOS resource-fork files and are not part of the analysis workflow.

## Main Workflow

The recommended execution order is:

1. Run `01.gene_selection_activate.ipynb` to generate WM-related gene scores and selected expression matrices.
2. Run `02.gene_coexprssion.ipynb` to construct transcriptional similarity matrices.
3. Run `03.corr.ipynb` to test the association between transcriptional similarity and F-TRACT response probability.
4. Run `03-1.control.ipynb` to perform spatial autocorrelation and random gene control analyses.
5. Run `04.virtual_knockout.ipynb` to calculate GCI values and split genes into GCI+ and GCI- groups.
6. Run `05.new_gene_score.ipynb` to prepare gene-score files for the GCI groups.
7. Run `06.cell_type_deconvolution.ipynb` and `10.enrichment.R` for cell-type and pathway enrichment analyses.
8. Run `07.PLS_posi_nega.ipynb` to compare GCI gene-expression profiles with receptor and macroscale brain maps.
9. Run `08.brainspan.ipynb` and `09.PCA.ipynb` to analyze developmental expression patterns.

## Data Sources

This workflow uses publicly available datasets and resources:

- Allen Human Brain Atlas (AHBA): http://human.brain-map.org
- Updated AHBA MNI coordinates: https://github.com/gdevenyi/AllenHumanGeneMNI
- Neurosynth working-memory activation map: https://neurosynth.org
- F-TRACT project: https://f-tract.eu/
- F-TRACT atlas release: https://zenodo.org/records/7015415
- neuromaps toolbox and reference brain maps: https://netneurolab.github.io/neuromaps/
- BrainSpan Atlas of the Developing Human Brain: https://www.brainspan.org
- Lausanne atlas resources described by Cammoun et al. (2012)

Raw datasets are not redistributed in this repository. Users should download data from the original sources and update local paths before running the scripts.

## Software Requirements

The analyses were performed using Python, R, and Jupyter notebooks. Major Python dependencies include:

- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `matplotlib`
- `plotnine`
- `nibabel`
- `nilearn`
- `brainsmash`
- `neuromaps`
- `pyls`
- `statsmodels`

The R enrichment analysis requires:

- `clusterProfiler`
- `ReactomePA`
- `org.Hs.eg.db`

Additional tools used in the broader workflow include `abagen`, `easy_lausanne`, and ENIGMA spin-test utilities.

## Path Configuration

Many scripts currently contain local absolute paths from the original analysis environment, for example paths beginning with:

```text
/Volumes/KWX/Grad/Programs/2025/202505/20250529/activity_gene_analysis/
```

Before reproducing the analysis, update these paths to match your local project directory. A recommended structure is:

```text
project_root/
  rawdata/
  processed/
  res/
  figures/
  code/
```

## Outputs

The workflow generates intermediate and final outputs including:

- WM-related gene-score files
- selected gene-expression matrices
- transcriptional similarity matrices
- aligned F-TRACT response probability matrices
- spatial and random-gene null distributions
- GCI+ and GCI- gene lists
- enrichment tables
- PLS results for receptor and macroscale brain maps
- BrainSpan developmental expression summaries
- figure source files and publication figures

## Notes for Reproducibility

- The F-TRACT response probability matrix is asymmetric; analyses should preserve directed regional pairs where specified.
- Negative WM activation values are set to zero before gene scoring.
- The main analysis uses the top 10% WM-matched genes and the top 15% left-hemisphere WM-activated regions.
- Spatial autocorrelation controls use BrainSMASH surrogate maps.
- PLS significance is assessed using spatially constrained permutation testing.
- Developmental expression analyses use four age periods: 2-8, 11-18, 19-30, and 36-40 years.

## Citation

If using this code, please cite the associated manuscript:

**Correlated Gene Expression Supports Signal Propagation in Human Working Memory Networks**.

## License

No license has been specified yet. Add a license file before public reuse if the repository is intended for open distribution.
