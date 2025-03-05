# mEdit

<!-- Badges -->
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/Interventional-Genomics-Unit/mEdit/main?logo=github)
![GitHub top language](https://img.shields.io/github/languages/top/Interventional-Genomics-Unit/mEdit?logo=github)
![PyPI - Version](https://img.shields.io/pypi/v/meditability)

<!-- Table of Contents -->
# Table of Contents

- [What is mEdit?](#what-is-medit)
  * [Program Structure](#program-structure)
  * [Features](#features)
- [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
  * [Running Tests](#running-tests)
- [Usage](#usage)
- [FAQ](#faq)
- [License](#license)
- [Contact](#contact)

## What is mEdit?
### Program Structure
<div align="center"> 
  <img src="src/infographics/new_medit_concept.png" alt="screenshot" />
</div>

### Features
 * Reference Human Genome
   * mEdit uses the RefSeq human genome reference GRCh38.p14
   * Alternatively, the user can provide a custom human assembly. [See [db_set](#1-database-setup) for details]
 * Alternative Genomes
   * mEdit can work with alternative genomes which are compared to the reference assembly
   * Pangenomes made public by the HPRC are built into mEdit and can be included in the analysis in 'standard' mode
 * Flexible editing tool selection
   * Several endonucleases and base-editors are built into mEdit and can be requested in any combination. [See options in [guide\_prediction](#3-guide-prediction)].
   * Custom editing tools can also be ingested by mEdit. [See how to format custom editors in [guide\_prediction](#3-guide-prediction)]

## Getting Started
### Prerequisites
 * The current version has 3 prerequisites:
   * [PIP](#pip)
   * [Anaconda](#anaconda)
   * [Mamba](#mamba)

#### PIP
  - Make sure `gcc` is installed
    ```
    sudo apt install gcc
    ```
  - Also make sure your pip up to date
    ```
    python -m pip install --upgrade pip
    ```
    * or: 
    ```
    apt install python3-pip
    ```

#### Anaconda
mEdit utilizes Anaconda to build its own environments under the hood. 
In the example below, we assume a Linux x86_64 system. 
For other sytems, follow the instructions on [this page](https://www.anaconda.com/docs/getting-started/miniconda/install).

- Install Miniconda:
  - Download the installer: 
    ```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash ~/Miniconda3-latest-Linux-x86_64.sh
    ```

  - Set up channel priority and update conda:  
    ```
    conda update --all
    conda config --set channel_priority strict
    ```

#### Mamba
  - The officially supported way of installing Mamba is through Miniforge.
  - The Miniforge repository holds the minimal installers for Conda and Mamba specific to conda-forge.
  - Example install
      ```
      wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
      bash Miniforge3-$(uname)-$(uname -m).sh
      ```
  - Important warning:
    - The supported way of using Mamba requires that no other packages are installed on the `base` conda environment
  - Additional information on how to operate Mamba:
    - [Mamba official page](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
    - [Miniforge repository](https://github.com/conda-forge/miniforge)

### Installation
 * mEdit is compatible with UNIX-based systems running on Intel processors and it's conveniently available via pyPI:

    ```
    pip install meditability
    ```

### Running Tests

 - As a Snakemake-based application, mEdit supports dry runs.
 - A dry run evaluates the presence of supporting data, and I/O necessary for each process without actually processing the run. 
 - All mEdit programs can be used called with the `--dry` option

## Usage

* To obtain information on how to run mEdit and view its programs, simply execute with the  `—-help` flag
    ```
    medit —-help
    ```
   
| Command            | Description                                                                                         |
|--------------------|-----------------------------------------------------------------------------------------------------|
| `db_set`          | Setup the necessary background data to run mEdit.                                                   |
| `list`            | Prints the current set of editors available on mEdit.                                              |
| `guide_prediction`| The core mEdit program finds potential guides for variants specified on the input by searching a diverse set of editors. |
| `offtarget`       | Predict off-target effects for the guides found.                                                    |



* There are four programs available in the current version
  * [db_set](#1-database-setup): Set up the necessary background data to run mEdit. This downloads ~7GB of data. 
  * [list](#2-editor-list): Prints the current set of editors available on mEdit.  
  * [guide\_prediction](#3-guide-prediction): This program scans for potential guides for variants specified on the input by searching a diverse set of editors.  
  * [offtarget](#4-off-target-analysis): Predicts off-target effect for the guides found  


### 1. **Database Setup**
 * Database Setup is used to retrieve the required information and datasets to run medit. The contents include the reference human genome, HPRC pangenome vcf files, Refseq, MANE, clinvar and more. See the database structure below.

    ```
    mEdit db_set [-h] [-d DB_PATH] [-l] [-c CUSTOM_REFERENCE] [-t THREADS]
    ```

#### **Parameters:**

#### Reference Database Pre-Processing
| Argument             | Description |
|----------------------|-------------|
| `-d DB_PATH`        | Path where the `mEdit_database` directory will be created ahead of the analysis. Requires ~7.5GB of in-disk storage. **[default: `./mEdit_database`]** |
| `-c CUSTOM_REFERENCE` | Path to a custom human reference genome in FASTA format. **Chromosome annotation must follow a `>chrN` format (case sensitive).** |
| `-t THREADS`        | Number of cores to use for parallel decompression of mEdit databases. |

### 2. **Editor List**
- In the current version there are 24 endonuclease editors and 29 base editor stored within medit. list  prints out a list of both base editors and endonuclease editors with the parameters used for guide prediction.
    ```
    mEdit list [-h] [-d DB_PATH]
    ```

#### **Parameters:**

#### Available Editors and Base Editors (BEs)

| Argument  | Description |
|-----------|-------------|
| `-d DB_PATH` | Path to the `mEdit_database` directory created using the `db_set` program. **[default: `./mEdit_database`]** |

**Output**;

```
Available endonuclease editors:  
-----------------------------  
name: spCas9  
pam, pam_is_first: NGG, False  
guide_len: 20  
dsb_position: -3  
notes: requirements work for SpCas9-HF1, eSpCas9 1.1,spyCas9  
5'-xxxxxxxxxxxxxxxxxxxxNGG-3'  
-----------------------------
```

### 3. **Guide Prediction**

- `guide_prediction` is the main program to search for guides given a list of variants. The pathogenic variants can be searched either from the ClinVar database or a _de novo_ variant (these must be provided as genomic coordinates. See `--qtype` option). 
- mEdit first generates variant incorporated gRNAs using the reference human genome. If the user chooses "fast" the search will end with the human reference genome. However if the user chooses “standard” or “vcf” the medit program will also go on to predict the impact of alternative genomic variants on either the pangenome or user provided vcf file.

```
mEdit guide_prediction [-h] -i QUERY_INPUT [-o OUTPUT] [-d DB_PATH] [-j JOBTAG] [-m {fast,standard,vcf}] [-v CUSTOM_VCF] [--qtype {hgvs,coord}] [--editor EDITOR_REQUEST]
                              [--be BE_REQUEST] [--cutdist CUTDIST] [--dry] [--pam PAM] [--guidelen GUIDE_LENGTH] [--pamisfirst] [--dsb_pos DSB_POSITION]
                              [--edit_win EDITING_WINDOW] [--target_base {A,C,G,T}] [--result_base {A,C,G,T}] [--cluster] [-p PARALLEL_PROCESSES] [--ncores NCORES]
                              [--maxtime MAXTIME]
```

#### **Parameters:**

#### **Input/Output Options**
| Argument          | Description |
|------------------|-------------|
| `-i QUERY_INPUT`  | Path to a plain text file containing the query (or set of queries) for mEdit analysis. See `--qtype` for formatting options. |
| `-o OUTPUT`      | Path to root directory where mEdit outputs will be stored. **[default: `mEdit_analysis_<jobtag>/`]** |
| `-d DB_PATH`     | Path to the `mEdit_database` directory created using the `db_set` program. **[default: `./mEdit_database`]** |
| `-j JOBTAG`      | Tag associated with the current mEdit job. A random jobtag is generated by default. |
---
#### **mEdit Core Parameters**
| Argument             | Description |
|---------------------|-------------|
| `-m {fast,standard,vcf}` | Mode option determining how mEdit runs: <br> **fast** - Uses one reference genome. <br> **standard** - Uses a reference genome and pangenomes. <br> **vcf** - Requires a custom VCF file. **[default: `standard`]** |
| `-v CUSTOM_VCF`     | Path to a gunzip-compressed VCF file for `vcf` mode. |
| `--qtype {hgvs,coord}` | Query type: <br> **hgvs** - Uses RefSeq ID + HGVS nomenclature. <br> **coord** - Uses hg38 1-based coordinates. **[default: `hgvs`]** |
| `--editor EDITOR_REQUEST` | Specifies the set of editors: <br> **clinical** - Uses clinically relevant editors. <br> **custom** - Requires `--pam`, `--pamisfirst`, `--guidelen`, `--dsb_pos`. |
| `--be BE_REQUEST` | Enables base editors: <br> **off** - Disables base editor search. <br> **default** - Uses ABE & CBE with `NGG` PAM and 4-8bp editing window. <br> **custom** - Requires `--pam`, `--guidelen`, `--edit_win`, `--target_base`, `--result_base`. |
| `--cutdist CUTDIST` | Maximum variant start position distance from the editor cut site. (Not available for base editors). **[default: `7`]** |
| `--dry` | Perform a dry run of mEdit. |
---
#### **Custom Editor Options**
| Argument        | Description |
|----------------|-------------|
| `--pam PAM`    | Specifies the PAM sequence for custom guide or base editor searches. |
| `--guidelen GUIDE_LENGTH` | Guide sequence length for custom endonuclease/base editor searches. |
| `--pamisfirst` | Indicates if the PAM is before the guide sequence. |
| `--dsb_pos DSB_POSITION` | Double-strand cut site relative to PAM. Example: `-3` for spCas9, `18,22` for Cas12. |
| `--edit_win EDITING_WINDOW` | Specifies editing window size (two comma-separated integers). Example: `"4,8"` for CBE. |
| `--target_base {A,C,G,T}` | Specifies the target base for base editor modification (e.g., `"A"` for ABE). |
| `--result_base {A,C,G,T}` | Specifies the base that the target base will be converted to (e.g., `"G"` for ABE). |
---
#### **SLURM Options**
| Argument               | Description |
|-----------------------|-------------|
| `--cluster`          | Request job submission through SLURM. **[default: `None`]** |
| `-p PARALLEL_PROCESSES` | Number of parallel processes for SLURM or local machine parallelization. **[default: `1`]** |
| `--ncores NCORES`    | Number of cores for each parallel process. **[default: `2`]** |
| `--maxtime MAXTIME`  | Maximum allowed time per parallel job. Format: `H:MM:SS`. Example: `"2:00:00"` for 2 hours. **[default: `1:00:00`]** |


### 4. **Off-target Prediction**

- The `offtarget` program applies [Guidescan2](https://github.com/pritykinlab/guidescan-cli) on the guides found in [guide\_prediction](#3-guide-prediction) and reports a summarized data set including the CFD score among other metrics.

    ```
    mEdit offtarget [-h] [--dry] [-o OUTPUT] [-d DB_PATH] -j JOBTAG [--select_editors SELECT_EDITORS] [--dna_bulge DNA_BULGE] [--rna_bulge RNA_BULGE]
                           [--max_mismatch MAX_MISMATCH] [--cluster] [-p PARALLEL_PROCESSES] [--ncores NCORES] [--maxtime MAXTIME]
    ```

#### **Parameters:**

#### **Input/Output Options**
| Argument               | Description |
|------------------------|-------------|
| `-o OUTPUT`           | Path to the root directory where `mEdit guide_prediction` outputs were stored. `"mEdit offtarget"` cannot operate if this path is incorrect. **[default: `mEdit_analysis_<jobtag>/`]** |
| `-d DB_PATH`          | Path to the `mEdit_database` directory created using the `db_set` program. **[default: `./mEdit_database`]** |
| `-j JOBTAG`           | Tag associated with the `"mEdit guide_prediction"` job. `"mEdit offtarget"` will use the `OUTPUT` option to access this `JOBTAG`. |
| `--select_editors SELECT_EDITORS` | Comma-separated list of editors to be analyzed for off-target effects. **[default: `all`]** |
| `--dna_bulge DNA_BULGE` | Sets the number of insertions in the off-target sequence. **[default: `0`]** |
| `--rna_bulge RNA_BULGE` | Sets the number of deletions in the off-target sequence. **[default: `0`]** |
| `--max_mismatch MAX_MISMATCH` | Maximum allowable number of mismatches in off-target analysis. **[default: `3`]** |
---
#### **SLURM Options**
| Argument              | Description |
|----------------------|-------------|
| `--cluster`         | Request job submission through SLURM. **[default: `None`]** |
| `-p PARALLEL_PROCESSES` | Number of parallel processes for SLURM or local machine parallelization. **[default: `1`]** |
| `--ncores NCORES`   | Number of cores for each parallel process. **[default: `2`]** |
| `--maxtime MAXTIME` | Maximum allowed time per parallel job. Format: `H:MM:SS`. Example: `"2:00:00"` for 2 hours. **[default: `1:00:00`]** |



## License
Copyright ©20xx [see Other Notes, below]. The Regents of the University of California (Regents). All Rights Reserved. Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 408, Berkeley, CA 94704-1362,  otl@berkeley.edu, for commercial licensing opportunities.

[Optional: Created by John Smith and Mary Doe, Department of Statistics, University of California, Berkeley.]

IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

## FAQ

## Cite us

## Contact
