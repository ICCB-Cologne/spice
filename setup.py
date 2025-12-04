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
        'numpy',
        'scipy',
        'seaborn',
        'tqdm',
        'fire',
        'CNSistent',
        'pandas',
        'pyyaml',
        'joblib',
        'ortools==9.8.3296',
        # needs to be installed using conda
        # 'openfst-python==1.8.2', 
    ],
    python_requires='>=3.8',
)
