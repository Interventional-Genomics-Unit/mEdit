# == Native Modules

# == Installed Modules
from setuptools import setup, find_packages
# == Project Modules


setup(
    name='meditability',
    version='0.1.3.9',
    description='',
    author='Interventional Genomics Unit',
    author_email='',
    entry_points={
        "console_scripts": [
            "medit = prog:main"
        ]
    },
    packages=find_packages("src"),
    package_dir={"": "src"},
    package_data={"smk.config": ["*.yaml"],
                  "smk.envs": ["*.yaml"],
                  "smk.pipelines": ["*.smk"]},
    include_package_data=True,
    install_requires=[
        'snakemake>=7.32.4',
        'wrapt>=1.15.0',
        'biopython>=1.81',
        'pyyaml>=6.0',
        'importlib-resources>=6.1.1',
        'pytz>=2023.3',
        'boto3>=1.28.57'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
    ],
)
