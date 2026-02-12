from setuptools import setup, find_packages

setup(
    name='spice',
    version='0.1',
    author='Tom L Kaufmann',
    description='SPICE: Selection Patterns In somatic Copy-number Events',
    packages=find_packages(),
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'spice=spice.cli:main',
        ],
    },
    install_requires=[
        # 'numpy==1.26.4',
        # 'pandas==2.2.3',
        'numpy',
        'pandas',
        'scipy',
        'seaborn',
        'tqdm',
        'fire',
        'pyyaml',
        'joblib',
        'ortools==9.8.3296',
    ],
    extras_require={
        'snakemake': ['snakemake>=7.0'],
        'preprocessing': ['CNSistent'],
    },
    python_requires='>=3.8',
)
