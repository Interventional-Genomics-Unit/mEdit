# mEdit

## Getting started
### Dependencies
- PIP vX.XX
  - Make sure `gcc` is installed
    - `sudo apt install gcc`

  - Keep your pip up to date
    - `python -m pip install --upgrade pip`
    - or: `apt install python3-pip`

- Anaconda vX.XX
  - Install Miniconda:
    - Download the installer: https://docs.conda.io/projects/miniconda/en/latest/
    - `bash bash Miniconda3-latest-<your-OS>.sh`
  - Set up your conda environment:
    - `conda update --all`
    - `conda config --set channel_priority strict`
- Mamba vX.XX
  - The officially supported way of installing Mamba is through Miniforge.
    - Important:
      - The supported way of using Mamba requires that no other packages are installed on the `base` conda environment
    - More information on: https://github.com/conda-forge/miniforge
      - `wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge-pypy3-<your-OS>.sh`
      - `bash Miniforge-pypy3-<your-OS>.sh`
- AWS CLI vX.XX
  - Make sure you are signed in with your AWS credentials:
  - More information: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
    ```
    curl "https://awscli.amazonaws.com/awscli-exe-<your-OS>.zip" -o "awscliv2.zip"
    unzip awscliv2.zip
    sudo ./aws/install 
    ```
    - 