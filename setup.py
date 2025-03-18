# == Native Modules
# read the contents of your README file
from pathlib import Path
# == Installed Modules
from setuptools import setup, find_packages
# == Project Modules
# == Readme setup
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='meditability',
    version='0.5.14',
    description='',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Daniel Bellieny, Taylor Hudson',
    author_email='bellieny@berkeley.edu, thudson@berkeley.edu',
    entry_points={
        "console_scripts": [
            "medit = prog:main",
            "ncbi_cross_db = prog:cross_db"
        ]
    },
    packages=find_packages("src"),
    package_dir={"": "src"},
    package_data={"smk.config": ["*.yaml"],
                  "smk.envs": ["*.yaml"],
                  "smk.pipelines": ["*.smk"]},
    include_package_data=True,
    install_requires=[
        'snakemake~=7.32.4',
        'pulp<2.8.0',
        'pandas>=2.1.1',
        'numpy~=1.26.0',
        'wrapt>=1.15.0',
        'biopython>=1.81',
        'pyyaml>=6.0',
        'importlib-resources>=6.1.1',
        'pytz>=2023.3',
        'boto3>=1.28.57',
        'requests>=2.32.3',
        'alive-progress>=3.1.5',
        'xmltodict>=0.13.0'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
    ],
)
