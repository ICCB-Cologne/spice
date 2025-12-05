# SPICE: Selection Patterns In somatic Copy-number Events

![](doc/logo_banner.png)

**SPICE**, Selection Patterns In somatic Copy-number Events, is an event-level framework that infers discrete copy-number events from allele-specific profiles.


## Installation


### Prerequisites
- Python >= 3.8
- medicc2 (including openfst)

### Install from pip/conda (recommended)
Coming soon!

### Install from source

1. Install MEDICC2 using conda/mamba:
```bash
conda install -c bioconda -c conda-forge medicc2
```

Or better directly create a new conda environment with MEDICC2 inside of it 
```bash
conda create -n spice_env -c conda-forge -c bioconda medicc2
conda activate spice_env
```

2. Clone the repository:
```bash
git clone https://tomkau@bitbucket.org/schwarzlab/spice.git
cd spice
```

3. Install in development mode:
```bash
pip install -e .
```

This will install SPICE and all its dependencies, and make the `spice` command available in your shell.

### Dependencies

SPICE automatically installs the following dependencies:
- numpy
- scipy
- pandas
- pyyaml
- joblib
- ortools (version 9.8.3296)


## Configuration

SPICE uses a configuration file for each run which are specified using the `--config` flag.
This means you can keep multiple configs (e.g., in `configs/`) and select them at runtime.

Parameters and directories not specified in the provided config file are taken from the default config file `default_config.yaml`.

### Minimal `config.yaml` override example
Each config must contain a name and the location of the input copy-number file like so:
```yaml
name: example_run
input_files:
   copynumber: data/example_data.tsv
```

For other parameters that can be modified, see `default_config.yaml`.

### Relative vs absolute paths

- `directories.*` entries (e.g., `data_dir`, `results_dir`, `log_dir`) as well as input files can be given as relative or absolute paths.
   - If relative, SPICE resolves them against the package root (i.e., the repository directory when installed in editable mode).
   - If absolute, SPICE uses them as-is.

### WGD Detection

SPICE supports two ways to determine WGD (whole genome duplication) status per sample. The pipeline branches on WGD status and uses different FSTs and neutral CN values accordingly.

- Provided status via `wgd_status` file:
   - Set `input_files.wgd_status` in your config to a TSV file.
   - The file must have two columns: first column is the sample identifier (used as index), second column named `wgd` with boolean values (`True`/`False`).
   - Example:
      ```tsv
      sample_id	wgd
      SA123	True
      SA456	False
      ```


- Inferred WGD status:
   - If `input_files.wgd_status` is missing or empty, SPICE infers WGD using copy-number data and the method specified by `params.wgd_inference_method`.
   - Supported values:
      - `major_cn`: heuristic whether at least half of the major copy-number is greater or equal to 2
      - `ploidy_loh`: PCAWG-style rule combining ploidy and LOH fraction

Notes
- WGD status impacts neutral CN values and constraint solving throughout the pipeline, so ensure this is set or inferred correctly.
- For haplotype-specific data, neutral CN is 1 (noWGD) vs 2 (WGD); for total CN, 2 vs 4 respectively.

### Sex (XY/XX) Detection

SPICE supports resolving sample sex (XY vs XX) either via a provided file or automatic inference. This affects handling of `chrX` and `chrY` in preprocessing and splitting.

- Provided status via `xy_samples` file:
   - Set `input_files.xy_samples` in your config to a TSV file.
   - The file must have two columns: first column is the sample identifier (used as index), second column named `xy` with boolean values (`True`/`False`) indicating XY (male) vs XX (female).
   - Example:
      ```tsv
      sample_id	xy
      SA123	True
      SA456	False
      ```

- Inferred XY status:
   - If `input_files.xy_samples` is missing or empty, SPICE infers XY by checking if any segments exist on chromosome `chrY` for a sample.

Effects
- For XY samples with haplotype-specific CN, the minor copy number of `chrX` and `chrY` is set to 0 during preprocessing and splitting.
- For XX samples, `chrY` is excluded (no segments on `chrY`).


## Usage

### Basic usage

```bash
spice --config <path/to/config>
```

For example

```bash
spice --config configs/example_config.yaml
```


Note that spice automatically deletes previous runs of the same name when it is rerun.

### General logic

spice runs 5 steps:
- `split`: Split haplotypes and preprocess input
- `all_solutions`: Enumerate all valid evolutionary paths
- `disambiguate`: Select best path using k-nearest neighbors
- `large_chroms`: Use MCMC sampling for chromosomes with many events
- `combine`: Combine all events into the final output

Additionally, there is an optional `plot` step to generate figures. It is not part of `all` and must be run explicitly.

For each step, nonWGD and WGD samples are treated separately and samples are split by chromosome and allele to give the file IDs "sample:chrom:allele". For each step, each ID is calculated separately and stored as separate files.

Intermediate files can be removed using

```bash
spice --clean --config <path/to/config>
```


### Preprocessing Step Details

The preprocessing step runs automatically (unless `--skip-preprocessing` is provided) and prepares the input for robust event inference. It performs:

- Data normalization: ensures chromosome names use `chr` prefix; converts starts/ends to integers and adjusts starts to 0-based.
- CN capping and filtering: caps copy numbers at 8; removes segments shorter than 1kb.
- WGD resolution: loads from `wgd_status.tsv` or infers as described in WGD Detection.
- Sex resolution: loads from `xy_samples.tsv` or infers by presence of `chrY`; for XY samples with haplotype-specific CN, sets minor CN of `chrX` and `chrY` to 0.
- Neighbor merging: merges adjacent segments with identical CNs to reduce fragmentation.
- Telomeres and centromeres: fills telomeric regions and optionally bins/unifies centromeres (can be skipped with `--pre-skip-centromeres`).
- MEDICC2 phasing: optional phasing of haplotypes; can be skipped with `--pre-skip-phasing`.
- Short arms and bounds: handles short arms and aligns segment ends to reference chromosome lengths.

Run control:
- Use `--skip-preprocessing` to bypass this step and proceed directly to `split`.


### Plotting

```bash
# Plot inferred events (optional)
# Per sample
spice plot --config <path/to/config> --plot-sample <SAMPLE_ID>
spice plot --config <path/to/config> --plot-sample <SAMPLE_ID> --plot-unit-size

# Per ID (format: sample:chr:cn_a|cn_b)
spice plot --config <path/to/config> --plot-id <sample:chr:allele>
```

Notes for plotting:
- The `plot` step requires `final_events.tsv` (produced by `combine`).
- Output PNGs are saved to `plot_dir/{name}/` (see `directories.plot_dir` in config; defaults to `plots/`).
- `--plot-unit-size` switches per-sample plots to unit-size segments.

For interactive exploration, see `notebooks/plotting.ipynb`.


### Run specific pipeline steps (Advanced)

You can run individual steps or combinations by specifying them as arguments:

```bash
# Run only input splitting
spice split --config <path/to/config>

# Run path enumeration and kNN selection
spice all_solutions disambiguate --config <path/to/config>

# Run MCMC sampling for large chromosomes
spice large_chroms --config <path/to/config>
```

### Parallel processing

Use multiple cores:
```bash
# Use 8 cores
spice --cores 8
```
While using multiple cores can technically make execution faster (especially in the case when spice takes a long time for single runs), it can also slow down execution when there are many entries to loop over.
We usually recommend to only use multiple cores for the `large_chroms` pipeline step (see above) as it takes the longest per sample.

Note that parallel processing will disable logging for the different subprocesses.

### Logging output

Control where logging output is sent with the `--log` flag:

* `--log terminal` (default): Writes logs to terminal only
* `--log file`: Writes logs to file only
* `--log both`: Writes logs to both terminal and file

When using `--log file` or `--log both`, logs are saved to the configured log directory from the config with a filename pattern: `{name}_{timestamp}.log`


## Input Format

SPICE expects tab-separated input files with copy-number segments. See example file `data/example_data.tsv`.

Required columns:
- `sample`: Sample identifier
- `chrom`: Chromosome name
- `start`: Segment start position
- `end`: Segment end position
- `cn_a`: Copy number for allele A (haplotype-specific) or total CN
- `cn_b`: Copy number for allele B (haplotype-specific only)

## Output

Results are saved in `results/{name}/`

The two main outputs are `summary.tsv` which is a summary over each ID (i.e. combination of sample, chromosome and allele) and the `final_events_df.tsv` which is a list of all inferred events.

### Intermediate files

The intermediate files are saved with separate directories for WGD and non-WGD profiles:
- `chrom_data_full/`: Preprocessed chromosome data
- `full_paths_single_solution/`: Chromosomes with unique solutions
- `full_paths_multiple_solutions/`: Chromosomes requiring kNN selection
- `knn_solved_chroms/`: Results from kNN selection
- `mcmc_solved_chroms_large/`: Results from MCMC sampling

Intermediate files can be removed using

```bash
spice --clean --config <path/to/config>
```

## Advanced Usage

### Using with Snakemake

For cluster execution, use the provided Snakefile:
```bash
snakemake --rerun-triggers mtime --slurm -j250 --keep-going
```

### Python API

You can also import and use SPICE functions directly in Python. Note that it is important to run `spice.load_config(config_file)` before any other spice imports
```python
import spice
spice.load_config(config_file)
from spice.pipeline import full_paths_from_graph_with_sv_wrapper
from spice.knn_graph import solve_with_knn
from spice.mcmc_for_large_chroms import mcmc_event_selection
```

## Citation

If you use SPICE in your research, please cite: [TODO]

## License

[TODO]

## Contact

For questions and issues, please contact: tom.kaufmann@iccb-cologne.org

