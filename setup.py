# == Native Modules

# == Installed Modules
from setuptools import setup, find_packages
# == Project Modules


setup(
    name='meditability',
    version='0.1.2',
    description='',
    author='Interventional Genomics Unit',
    author_email='',
    entry_points={
        "console_scripts": [
            "medit = prog:main"
        ]
    },
    packages=find_packages(),
    install_requires=[
        'snakemake>=7.32.4',
        'biopython>=1.81',
        'pyyaml>=6.0',
        'pytz>=2023.3'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
    ],
)
