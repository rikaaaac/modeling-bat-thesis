## Folder structure (on gDrive)

- Thesis
  - README.md
  - deliverables
  - scrna-seq
    - README.md
    - codes
    - all h5ad files
    - results
      - regular (master + pdgfra subsets): all plots, all degs
      - trajectory: scvelo, cellrank, combined modalities
  - spatial
    - README.md
    - codes
    - all h5ad files
    - results: sec_cells, pdgfra, tdtomato

### How to navigate the folders

- All deliverables are in `deliverables/`. This contains SRI poster, my thesis final presentation slides, and my final thesis.
- All analyses relating to scRNA-seq are in `scrna-seq/`. It contains all codes I've written in one folder called `codes/`, all results together (with subfolders) in `results/`,
  and all intermediate and final h5ad files in `h5ad/`.
  - Read the `README.md` file within this folder for all explanations of codes, results, and various files. It also entailed the methods description.
- All analyses relating to visium/spatial transcriptomics are in `spatial`. It contains all codes I've written in one folder called `codes/`, all results together (with subfolders) in `results/`,
  and all intermediate and final h5ad files in `h5ad/`.
  - Read the `README.md` file within this folder for all explanations of codes, results, and various files. It also entailed the methods description.

## Folder structure (on github)

- scrna-seq
  - codes: all codes for scrna-seq analysis
  - README.md: detailed info for all code files
- spatial
  - codes: all codes for spatial transcriptomics analysis
  - README.md: detailed info for all code files
- thesis_figures_cleanup.ipynb: misc file when I was prettifying the figures for thesis write-up.

**Note:** for all results, deliverables, and other files, go to gDrive. Github repo is just for all codes.
