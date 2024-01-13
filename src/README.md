# mEdit

## Getting started
### Dependencies
- PIP
  - Make sure `gcc` is installed
    - `sudo apt-get install gcc`

  - Keep your pip up to date
    - `python -m pip install --upgrade pip`

- Anaconda
  - Set up your conda environment:
    - `conda update --all`
    - `conda config --set channel_priority strict`
- Mamba
  - Install mamba through conda 
    - `conda install -n base -c conda-forge mamba`
- AWS CLI
  - Make sure you are signed in with your AWS credentials
    ```
    curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
    unzip awscliv2.zip
    sudo ./aws/install 
    ```
    - 