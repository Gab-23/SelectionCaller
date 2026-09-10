# Alteration-type-agnostic detection of selection in cancer cohorts using machine learning

![SelectionCaller](SelectionCaller.png)

<!-- Badges - uncomment/edit as needed
[![CI](https://github.com/<org>/<repo>/actions/workflows/ci.yml/badge.svg)](https://github.com/<org>/<repo>/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-see%20LICENSE.md-blue)](LICENSE.md)
-->

## Table of contents

- [Overview](#overview)
- [Background](#background)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline modules](#pipeline-modules)
- [Testing and continuous integration](#testing-and-continuous-integration)
- [Running on a SLURM cluster](#running-on-a-slurm-cluster)
- [Outputs](#outputs)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Overview

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor
incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis
nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu
fugiat nulla pariatur.

## Background

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor
incididunt ut labore et dolore magna aliqua.

- Lorem ipsum dolor sit amet
- Consectetur adipiscing elit
- Sed do eiusmod tempor incididunt

## Repository structure

```
SelectionCaller/
├── .github/                # GitHub Actions CI workflows
├── .snakemake/             # Snakemake's own run metadata (not tracked)
├── .tests/                 # Unit tests and toy-data integration tests
│   ├── unit/
│   └── integration/
├── config/                 # Production configuration files (config.yaml, sample sheets, ...)
├── dags/                   # Lorem ipsum dolor sit amet
├── doc/                    # Lorem ipsum dolor sit amet
├── launcher/                # SLURM launch script(s)
├── resources/               # Intermediate, per-step data (gitignored)
├── results/                 # Plots and QC reports (gitignored)
├── temp_DARs/                # Lorem ipsum dolor sit amet
├── temp_ML/                  # Lorem ipsum dolor sit amet
├── temp_plots/                # Lorem ipsum dolor sit amet
├── workflow/
│   ├── envs/                   # Per-rule conda environment files
│   ├── rules/                  # Snakemake rule files, one module per subdirectory
│   ├── scripts/                 # R/Python analysis scripts called by the rules
│   ├── snakefile_helpers/        # Shared helper functions imported by the Snakefiles
│   ├── snakefile_bin_level        # Main Snakefile: bin-level genomic feature generation
│   ├── snakefile_ML               # Lorem ipsum dolor sit amet
│   └── snakefile_peak_level       # Lorem ipsum dolor sit amet
├── .gitignore
├── LICENSE.md
├── README.md
└── SelectionCaller.png
```

## Requirements

Lorem ipsum dolor sit amet, consectetur adipiscing elit:

- Snakemake (>= lorem ipsum version)
- Conda/Mamba
- Lorem ipsum dolor sit amet

## Installation

```bash
git clone <repository-url>
cd SelectionCaller
# Lorem ipsum dolor sit amet
```

## Configuration

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor
incididunt ut labore et dolore magna aliqua. See `config/config.yaml` for the
full list of parameters.

## Usage

```bash
snakemake --snakefile workflow/snakefile_bin_level \
          --configfile config/config.yaml \
          --cores <N> \
          --use-conda
```

Lorem ipsum dolor sit amet, consectetur adipiscing elit.

## Pipeline modules

The `snakefile_bin_level` pipeline generates a genome-wide, bin-level feature
matrix from several independent data modalities:

| Module | Description |
|---|---|
| Backbone/annotations | Lorem ipsum dolor sit amet |
| DepMap | Lorem ipsum dolor sit amet |
| IntOGen | Lorem ipsum dolor sit amet |
| Epigenomes (Roadmap) | Lorem ipsum dolor sit amet |
| COSMIC | Lorem ipsum dolor sit amet |
| GDC (copy number) | Lorem ipsum dolor sit amet |
| Fragile sites | Lorem ipsum dolor sit amet |
| Ohnologs/paralogs | Lorem ipsum dolor sit amet |
| Gene density | Lorem ipsum dolor sit amet |
| Haploinsufficiency | Lorem ipsum dolor sit amet |
| Reference genome | Lorem ipsum dolor sit amet |
| RepliSeq | Lorem ipsum dolor sit amet |
| ATAC-seq | Lorem ipsum dolor sit amet |
| Hi-C | Lorem ipsum dolor sit amet |
| Combined feature matrix | Lorem ipsum dolor sit amet |

`snakefile_ML` and `snakefile_peak_level` — lorem ipsum dolor sit amet,
consectetur adipiscing elit.

## Testing and continuous integration

This repository ships with two layers of automated testing, run on every push
and pull request via GitHub Actions (see `.github/workflows/ci.yml`):

- **Unit tests** (`.tests/unit/`) — exercise individual R/Python scripts against
  small, hand-built toy inputs, independently of Snakemake.
- **Toy integration test** (`.tests/integration/`) — runs the actual Snakemake
  pipeline end-to-end on a small, fake genome assembly with synthetic data, to
  catch wiring/DAG-level issues that unit tests alone cannot.

To run the tests locally:

```bash
# unit tests
for f in .tests/unit/*.sh; do bash "$f"; done

# dry-run against the production config
snakemake -n -s workflow/snakefile_bin_level --configfile config/config.yaml
```

Lorem ipsum dolor sit amet, consectetur adipiscing elit.

## Running on a SLURM cluster

A ready-to-use launcher script is provided under `launcher/`. Submit it with:

```bash
sbatch launcher/<script_name>.sh workflow/snakefile_bin_level config/config.yaml <cores>
```

Lorem ipsum dolor sit amet, consectetur adipiscing elit.

## Outputs

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor
incididunt ut labore et dolore magna aliqua.

## Citation

Lorem ipsum dolor sit amet, consectetur adipiscing elit.

```
Lorem ipsum citation placeholder.
```

## License

See [LICENSE.md](LICENSE.md).

## Contact

Lorem ipsum dolor sit amet, consectetur adipiscing elit.