# GPCRchimeraDB

Main repository of the database **GPCRchimeraDB** ([https://www.bio2byte.be/gpcrchimeradb/](https://www.bio2byte.be/gpcrchimeradb/)).

## Format the Data

In the folder `Notebooks`, the following Jupyter notebooks are available:
- [**Format a chimera's data**](Notebooks/create_chimera.ipynb): Identify parents, cutting sites, and mutations based on its sequence.
- [**Generate entry file for chimeric GPCRs**](Notebooks/json_GPCRchimeraDB_chimera.ipynb).
- [**Generate entry file for natural GPCRs**](Notebooks/json_GPCRchimeraDB_natural.ipynb).

## The Data

Available on Zedodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10854343.svg)](https://doi.org/10.5281/zenodo.14989075) :
- Entry files (JSON, 1 per entry).
- FASTA sequences.
- Structures (PDB, AF2, AF2 MS, ESMFold, dssp, mapping).
- Master alignment (MSA).
- Activation state structures (GPCRdb).
- Classification.
- G-protein and B-arrestin coupling (GPCRdb).
- Endogenous ligands (GPCRdb).
- Variants (GPCRdb).
- Scientific to common names organisms (UniProt).

## The Database

**Explore your favorite chimeric or natural GPCR at the sequence, structural, and biophysical levels and design new chimeras!**

Visit GPCRchimeraDB: [https://www.bio2byte.be/gpcrchimeradb/](https://www.bio2byte.be/gpcrchimeradb/).

## The Database's Documentation

Visit GPCRchimeraDB's documentation: [https://gpcrchimeradb-docs.readthedocs.io/en/latest/index.html](https://gpcrchimeradb-docs.readthedocs.io/en/latest/index.html).

## Citation

```bibtex
@article{CRAUWELS2025169164,
title = {GPCRchimeraDB: A database of chimeric G protein-coupled receptors (GPCRs) to assist their design},
journal = {Journal of Molecular Biology},
pages = {169164},
year = {2025},
issn = {0022-2836},
doi = {https://doi.org/10.1016/j.jmb.2025.169164},
author = {Charlotte Crauwels and Adrián Díaz and Wim Vranken},
}
```


## Contact Us

[bio2byte@vub.be](mailto:bio2byte@vub.be)

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
