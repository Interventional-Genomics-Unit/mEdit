# mEdit

<!-- Badges -->
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/Interventional-Genomics-Unit/mEdit/total?logo=github)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/Interventional-Genomics-Unit/mEdit/main?logo=github)
![GitHub top language](https://img.shields.io/github/languages/top/Interventional-Genomics-Unit/mEdit?logo=github)
![PyPI - Version](https://img.shields.io/pypi/v/meditability)
![GitHub License](https://img.shields.io/github/license/Interventional-Genomics-Unit/mEdit)

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
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [FAQ](#faq)
- [License](#license)
- [Contact](#contact)

## What is mEdit?
### Program Structure
<div align="center"> 
  <img src="src/infographics/mEdit_graphical_overview.png" alt="screenshot" />
</div>

### Features
 * Reference Human Genome
   * mEdit has three different options for loading reference genomes:
   * 1- The latest compatible reference used by the HPRC on human pangenome assemblies
   * 2- The latest RefSeq human genome reference
   * 3- Custom version provided by the user
## Getting Started
### Prerequisites
 * The current version has 4 prerequisites:
   * [PIP](#pip)
   * [Anaconda](#anaconda)
   * [Mamba](#mamba)
   * [AWS-CLI](#aws-cli)
   * [Tabix](#tabix)
   
#### PIP
  - Make sure `gcc` is installed
    ```
    sudo apt install gcc
    ```
  - Make sure your pip up to date
    ```
    python -m pip install --upgrade pip
    ```
    * or: 
    ```
    apt install python3-pip
    ```
#### Anaconda
  - Install Miniconda:
    * Download the installer at: https://docs.conda.io/projects/miniconda/en/latest/ 
    ```
    bash Miniconda3-latest-<your-OS>.sh
    ```
  - Set up your conda environment: 
    ```
    conda update --all
    conda config --set channel_priority strict
    ```
 
#### Mamba
  - The officially supported way of installing Mamba is through Miniforge.
    * Important:
      * The supported way of using Mamba requires that no other packages are installed on the `base` conda environment
    * More information about miniforge: https://github.com/conda-forge/miniforge
    * Details on how to correctly install Mamba: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
      ```
      wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge-pypy3-<your-OS>.sh
      bash Miniforge-pypy3-<your-OS>.sh
      ```
#### AWS CLI
  - Make sure you are signed in with your AWS credentials:
  - More information: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
    ```
    curl "https://awscli.amazonaws.com/awscli-exe-<your-OS>.zip" -o "awscliv2.zip"
    unzip awscliv2.zip
    sudo ./aws/install 
    ```
#### Tabix
 * This allows mEdit to use bgzip for file compression.
   * In linux systems, install `tabix`.
      ```
        apt install tabix
     ```
   * On macOS machines, intall `htslib`.
      ```
        brew install htslib
     ```

### Installation
mEdit is compatible with UNIX-based systems running on Intel processors and it's conveniently available via pyPI:
```
pip install meditability
```

### Running Tests

## Usage
  - `db_set`
  - `list`
  - `guide_prediction`
  - `offtargets`
## Roadmap
## Contributing
## FAQ
## License
## Contact
