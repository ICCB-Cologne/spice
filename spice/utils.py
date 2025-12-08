import logging
import os
import pickle
import re
import itertools
import warnings
from functools import wraps
from datetime import datetime

from tqdm.auto import tqdm
import pandas as pd
import numpy as np
import yaml

from scipy.stats import spearmanr, permutation_test

# Global logging configuration (set by CLI or defaults)
_LOGGING_CONFIG = {
    'mode': 'terminal',  # 'terminal', 'file', or 'both'
    'log_dir': None,
    'config_name': None,
    'level': 'INFO',
    'file_handler': None,  # Shared file handler for all loggers
    'configured': False
}

def configure_logging(log_mode, log_dir, config_name, level: str = 'INFO'):
    """
    Configure global logging settings for all SPICE loggers.
    Should be called once by CLI before any logging occurs.
    
    Args:
        log_mode: 'terminal', 'file', or 'both'
        log_dir: Directory for log files
        config_name: Name from config (used in log filename)
        steps: List of steps being run (used in log filename)
    """
    global _LOGGING_CONFIG
    
    _LOGGING_CONFIG['mode'] = log_mode
    _LOGGING_CONFIG['log_dir'] = log_dir
    _LOGGING_CONFIG['config_name'] = config_name
    _LOGGING_CONFIG['level'] = level
    _LOGGING_CONFIG['configured'] = True
    
    # Create shared file handler if needed
    if log_mode in ['file', 'both']:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_filename = f"{config_name}_{timestamp}.log"
        log_path = os.path.join(log_dir, log_filename)
        
        os.makedirs(log_dir, exist_ok=True)
        file_handler = logging.FileHandler(log_path)
        formatter = logging.Formatter('*** %(name)s - %(levelname)s - %(asctime)s - %(message)s')
        file_handler.setFormatter(formatter)
        _LOGGING_CONFIG['file_handler'] = file_handler
        
        if log_mode == 'file':
            print(f"Logging to file: {log_path}")


def get_logger(name, load_config=True, config_file=None):
    """
    Get a logger that respects global logging configuration.
    
    Args:
        name: Logger name
        load_config: Whether to load logging level from config (for backward compatibility)
        config_file: Config file to load (for backward compatibility)
    
    Returns:
        Configured logger instance
    """
    global _LOGGING_CONFIG
    
    # Determine logging level
    if load_config:
        if config_file is None:
            if os.path.exists('../config.yaml'):
                config_file = '../config.yaml'
            elif os.path.exists('../default_config.yaml'):
                config_file = '../default_config.yaml'
    
    if _LOGGING_CONFIG['configured'] and _LOGGING_CONFIG.get('level'):
        logging_level = _LOGGING_CONFIG['level']
    else:
        if config_file is not None:
            with open(config_file, 'rt') as f:
                config = yaml.safe_load(f.read())
            logging_level = config['params']['logging_level']
        else:
            logging_level = "INFO"
    
    logger = logging.getLogger(name)
    logger.setLevel(logging_level)
    logger.handlers.clear()  # Clear any existing handlers
    
    formatter = logging.Formatter('*** %(name)s - %(levelname)s - %(asctime)s - %(message)s')
    
    # If global config is set, use it; otherwise default to terminal
    log_mode = _LOGGING_CONFIG['mode'] if _LOGGING_CONFIG['configured'] else 'terminal'
    
    # Add terminal handler if needed
    if log_mode in ['terminal', 'both']:
        terminal_handler = logging.StreamHandler()
        terminal_handler.setFormatter(formatter)
        logger.addHandler(terminal_handler)
    
    # Add file handler if needed and available
    if log_mode in ['file', 'both'] and _LOGGING_CONFIG['file_handler'] is not None:
        logger.addHandler(_LOGGING_CONFIG['file_handler'])
    
    return logger


def resolve_data_file(return_raw=False) -> str:
    """Resolve the chromosome segments file path.
    """
    from spice import config, directories
    logger = get_logger('utils')

    name = config.get('name')
    data_dir = config['directories']['data_dir']
    orig = config['input_files']['copynumber']
    processed = os.path.join(data_dir, f"{name}_processed.tsv")
    if not return_raw:
        orig = orig.replace('.tsv', '_split.tsv')
        processed = os.path.join(data_dir, f"{name}_processed_split.tsv")

    # Prefer processed split file if available, else original
    cur_file = orig
    if processed and os.path.exists(processed):
        cur_file = processed
    
    if not os.path.isabs(cur_file):
        cur_file = os.path.join(directories['base_dir'], cur_file)
    log_debug(logger, f"Resolved chrom_segments_file: {cur_file}")
    return cur_file


def load_final_events():
    from spice import config, directories
    from spice.utils import get_logger
    logger = get_logger('utils')
    results_dir = os.path.join(directories['results_dir'], config['name'])
    if not os.path.exists(os.path.join(results_dir, 'final_events.tsv')):
        raise FileNotFoundError(f"final_events.tsv not found in dir {results_dir}. Run SPICE event inference first")
    final_events_df = pd.read_csv(
        os.path.join(results_dir, 'final_events.tsv'), sep='\t', dtype={'cn': str, 'diff': str})
    return final_events_df


def step_aware_cleanup(results_dir, requested_steps=None):
    """Delete artifacts from the first requested step onward, including preprocessing.

    Steps order: ['preprocessing', 'split', 'all_solutions', 'disambiguate', 'large_chroms', 'combine']
    - Preprocessing artifacts under `directories.data_dir`:
      - `{name}_processed.tsv`
      - `{name}_processed_split.tsv`
    - Other artifacts under `results_dir` per WGD branch.
    """
    import shutil
    from spice import config, directories
    from spice.utils import get_logger
    logger = get_logger('utils')
    name = config.get('name')

    step_order = ['preprocessing', 'split', 'all_solutions', 'disambiguate', 'large_chroms', 'combine']
    if requested_steps and all(s not in step_order for s in (requested_steps or [])):
        log_debug(logger, 'No step-aware cleanup requested (no --keep-old).')
        return

    step_targets = {
        'preprocessing': [
            ('__data__', f'{name}_processed.tsv'),
            ('__data__', f'{name}_processed_split.tsv'),
        ],
        'split': [
            ('nowgd', 'chrom_data_full'),
            ('nowgd', 'chrom_data_large'),
            ('wgd', 'chrom_data_full'),
            ('wgd', 'chrom_data_large'),
            ('', 'WGD_status_major_cn.png'),
            ('', 'WGD_status_ploidy_loh.png'),
        ],
        'all_solutions': [
            ('nowgd', 'full_paths_multiple_solutions'),
            ('nowgd', 'full_paths_single_solution'),
            ('wgd', 'full_paths_multiple_solutions'),
            ('wgd', 'full_paths_single_solution'),
        ],
        'disambiguate': [
            ('nowgd', 'knn_solved_chroms'),
            ('wgd', 'knn_solved_chroms'),
        ],
        'large_chroms': [
            ('nowgd', 'mcmc_solved_chroms_large'),
            ('wgd', 'mcmc_solved_chroms_large'),
        ],
        'combine': [
            ('', 'final_events.tsv'),
            ('', 'summary.tsv'),
            ('', 'failed_reports.tsv'),
        ],
    }

    def delete_path(rel_root: str, rel_name: str):
        if rel_root == '__data__':
            data_dir = config['directories']['data_dir']
            target = os.path.join(data_dir, rel_name) if data_dir else None
        else:
            target = os.path.join(results_dir, rel_root, rel_name) if rel_root else os.path.join(results_dir, rel_name)
        if not target:
            return
        if os.path.isdir(target):
            log_debug(logger, f'Removing directory {target}')
            shutil.rmtree(target, ignore_errors=True)
        elif os.path.isfile(target):
            log_debug(logger, f'Removing file {target}')
            try:
                os.remove(target)
            except FileNotFoundError:
                pass

    start_idx = 0 if requested_steps is None else min(step_order.index(s) for s in requested_steps if s in step_order)

    targets = []
    for s in step_order[start_idx:]:
        targets.extend(step_targets.get(s, []))

    if targets:
        logger.info('Step-aware cleanup enabled (no --keep-old).')
        for rel_root, rel_name in targets:
            delete_path(rel_root, rel_name)
    else:
        logger.info('No mapped artifacts to clean for requested steps.')


def chrom_id_from_id(cur_id):
    return re.sub(':cn_[ab](:\d+$)?', '', cur_id)


def get_sister_allele(cur_id):
    if "cn_a" in cur_id:
        return cur_id.replace("cn_a", "cn_b")
    else:
        return cur_id.replace("cn_b", "cn_a")

class CALC_NEW_SKIP():
    def __repr__(self):
        return 'CALC_NEW skip signal'
    def __eq__(self, other):
        return other is None


class CALC_NEW:
    def __init__(self, filename=None, force_new=False, verbose=True):
        self.filename = filename
        self.force_new = force_new
        self.logger = logging.getLogger('CALC_NEW')
        if verbose:
            self.logger.setLevel(logging.INFO)
        else:
            self.logger.setLevel(logging.WARNING)


    def __call__(self, func):
        def wrapper(*args, **kwargs):

            self.result = None
            self.save_data = True

            if 'calc_new_filename' in kwargs:
                filename = kwargs['calc_new_filename']
                del kwargs['calc_new_filename']
            else:
                filename = self.filename
            if 'calc_new_force_new' in kwargs:
                force_new = kwargs['calc_new_force_new']
                del kwargs['calc_new_force_new']
            else:
                force_new = self.force_new
            if 'calc_new_verbose' in kwargs:
                if kwargs['calc_new_verbose']:
                    self.logger.setLevel(logging.INFO)
                else:
                    self.logger.setLevel(logging.WARNING)
                del kwargs['calc_new_verbose']

            if filename is None:
                log_debug(self.logger, f"File path is None, ignoring calc_new")
                self.result = func(*args, **kwargs)
                return self.result
            else:
                if os.path.exists(filename):
                    if not force_new:
                        log_debug(self.logger, f"Data exists at '{filename}'")
                        self.save_data = False
                        log_debug(self.logger, f"Loading data")
                        with open(filename, 'rb') as f:
                            self.result = pickle.load(f)
                        return self.result
                    else:
                        log_debug(self.logger, f"Data exists at but forcing recalculation for '{filename}'")
                        self.result = None
                else:
                    log_debug(self.logger, f"Data does not exist at {filename}")

                if self.result is None:
                    self.result = func(*args, **kwargs)
                    if isinstance(self.result, CALC_NEW_SKIP):
                        log_debug(self.logger, f"Function returned skip signal, not saving data")
                        self.save_data = False
                        self.result = None
                    else:
                        log_debug(self.logger, f"Save data to {filename}")
                        if not os.path.exists(os.path.dirname(filename)):
                            os.makedirs(os.path.dirname(filename), exist_ok=True)
                        with open(filename, 'wb') as f:
                            pickle.dump(self.result, f)

                return self.result
        return wrapper


def assert_close(a, b):
    assert np.isclose(a, b), print(f"{a} - {b} = {a - b}")


def open_pickle(filename, n_elements=None, fail_if_nonexisting=True,
                fail_for_other_errors=True, data_type=None):

    if not fail_if_nonexisting and not os.path.exists(filename):
        local_logger = get_logger('utils')
        local_logger.warning('File does not exist, returning None')
        return None
    try:
        with open(filename, 'rb') as f:
            data = pickle.load(f)
    except Exception as e:
        if isinstance(e, FileNotFoundError):
            if fail_if_nonexisting:
                raise Exception(f"File {filename} was not found:\n{str(e)}")
        elif fail_for_other_errors:
            raise Exception(f"Opening pickle {filename} failed with error code:\n{str(e)}")
        data = None
    if n_elements is not None and not hasattr(data, '__len__') and data is None:
        data = tuple([None] * n_elements)

    if data_type is not None:
        assert isinstance(data, data_type), f"Data type is {type(data)}, expected {data_type}"

    return data


def save_pickle(obj, filename, create_dir_if_not_exists=True):

    if create_dir_if_not_exists and not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename), exist_ok=True)

    with open(filename, 'wb') as f:
        pickle.dump(obj, f)


def linkage_order(data):
    import scipy.cluster.hierarchy as shc
    linkage = shc.linkage(data, method='ward')
    dend = shc.dendrogram(linkage, no_plot=True)

    return dend['leaves']


# TODO: Create to events_from_graph?
def create_full_df_from_diff_df(df, cur_id, cur_chrom_segments, chrom_lengths=None,
                                calc_telomere_bound=True):
    if chrom_lengths is None:
        from spice.data_loaders import load_chrom_lengths
        chrom_lengths = load_chrom_lengths()
    cur_chrom_segments = cur_chrom_segments.query('id == @cur_id')

    df = df.copy()
    df['id'] = cur_id
    df[['sample', 'chrom', 'allele']] = cur_id.split(':')
    df['chrom_length'] = chrom_lengths.loc[df['chrom'].values[0]]
    df['chrom_id'] = cur_id[:cur_id.rfind(':')]
    df['start'] = cur_chrom_segments['start'].values[df['diff'].map(lambda x: x.find('1')).values]
    df['end'] = cur_chrom_segments['end'].values[df['diff'].map(lambda x: x.rfind('1')).values]
    # width cannot simply be end - start because of LOH-shortened events
    df['total_loh_gaps'] = [cur_chrom_segments.eval('end - start').values[
        diff.find('1'):diff.rfind('1')][
            np.fromiter(diff[diff.find('1'):diff.rfind('1')], int) == 0].sum() for diff in df['diff'].values]
    df['width'] = df.eval('(end - start) - total_loh_gaps')
    
    if calc_telomere_bound:
        df[['telomere_bound', 'whole_arm', 'whole_chrom']] = np.stack(calc_telomere_bound_whole_arm_whole_chrom(df), axis=1)
       
    if 'is_gain' in df.columns:
        df['type'] = df['is_gain'].map({True: 'gain', False: 'loss'})
    if 'type' in df.columns:
        df['type_diff'] = df[['type', 'diff']].apply(lambda x: f'{x["type"]}:{x["diff"]}', axis=1)

    return df


def calc_telomere_bound_whole_arm_whole_chrom(data, return_left_and_right=False):
    telomere_bound_l, telomere_bound_r = calc_telomere_bound_left_and_right(data)
    centromere_bound_l, centromere_bound_r = calc_centromere_bound(data)

    # Fix problems with chrX assignment
    if isinstance(data, pd.DataFrame):
        left_bool = data.eval('chrom == "chrX" and start < 28e5').values
        right_bool = data.eval('chrom == "chrX" and end > 1548e5').values
        telomere_bound_l = np.logical_or(telomere_bound_l, left_bool)
        telomere_bound_r = np.logical_or(telomere_bound_r, right_bool)
    # else:
    #     logger.warning("Cannot fix chrX telomere bounds for non-DataFrame data")
    telomere_bound = np.logical_or(telomere_bound_l, telomere_bound_r)
    whole_chrom = np.logical_and(telomere_bound_l, telomere_bound_r)
    whole_arm = np.logical_and(
        (~whole_chrom),
        np.logical_or(
            np.logical_and(telomere_bound_l, centromere_bound_r),
            np.logical_and(telomere_bound_r, centromere_bound_l)
        )
    )
    if return_left_and_right:
        return centromere_bound_l, centromere_bound_r, telomere_bound_l, telomere_bound_r, telomere_bound, whole_arm, whole_chrom
    else:
        return telomere_bound, whole_arm, whole_chrom


def calc_telomere_bound_left_and_right(data):
    if isinstance(data, pd.DataFrame):
        telomere_bound_l = data['diff'].str.startswith('1').values
        telomere_bound_r = data['diff'].str.endswith('1').values

    elif isinstance(data, tuple):
        assert len(data) == 6, f'Data tuple must have 6 elements. Current has {len(data)} elements.'
        event_sorted, cn_profile, starts, ends, centro_start, centro_end = data
        telomere_bound_l = (event_sorted == 0).any(axis=1)
        telomere_bound_r = (event_sorted == len(cn_profile)).any(axis=1)
    else:
        raise ValueError(f"Data type {type(data)} not supported")

    return telomere_bound_l, telomere_bound_r


def calc_centromere_bound(data):
    if isinstance(data, pd.DataFrame):
        from spice.data_loaders import load_centromeres # import has to be here to prevent circular imports
        centromeres = load_centromeres().astype(int)
        data = (data
                .drop(columns=['centro_start', 'centro_end'], errors='ignore')
                .join(centromeres, on='chrom'))

        centromere_bound_r = data.eval('end >= centro_start-2 and end <= centro_end+2').values
        centromere_bound_l = data.eval('start >= centro_start-2 and start <= centro_end+2').values

    elif isinstance(data, tuple):
        assert len(data) == 6, f'Data tuple must have 6 elements. Current has {len(data)} elements.'
        event_sorted, cn_profile, starts, ends, centro_start, centro_end = data
        # Note: safer to use np.logical_and and np.logical_or instead of & and |
        centromere_bound_r = np.logical_and(ends >= centro_start-2, ends <= centro_end+2)
        centromere_bound_l = np.logical_and(starts >= centro_start-2, starts <= centro_end+2)
    else:
        raise ValueError(f"Data type {type(data)} not supported")

    return centromere_bound_l, centromere_bound_r


@CALC_NEW()
def create_chrom_type_pos_indices(events_df):
    assert (events_df.index == np.arange(len(events_df))).all(), "index has to be 0, 1, 2, 3, ..."
    chrom_type_pos_indices = dict()
    for chrom, type, pos in itertools.product(events_df['chrom'].unique(), events_df['type'].unique(), events_df['pos_detail'].unique()):
        chrom_type_pos_indices[(chrom, type, pos)] = events_df.query('chrom == @chrom and type == @type and pos_detail == @pos').index.values

    assert len(np.unique(np.concatenate([x for x in chrom_type_pos_indices.values()]))) == len(events_df), f"{len(np.unique(np.concatenate([x for x in chrom_type_pos_indices.values()]))), len(events_df)}"
    assert all(np.sort(np.unique(np.concatenate([x for x in chrom_type_pos_indices.values()]))) == np.sort(events_df.index))

    return chrom_type_pos_indices


def create_full_paths_events_df(cur_full_paths, chrom_segments):
    unique_events_df = pd.DataFrame(cur_full_paths.events.values())
    unique_events_df = create_full_df_from_diff_df(unique_events_df, cur_full_paths.id, chrom_segments)
    cur_events_df = pd.concat([unique_events_df.iloc[np.array(list(sol.elements()))] for sol in cur_full_paths.solutions], axis=0)
    cur_events_df['chain_nr'] = np.repeat(np.arange(cur_full_paths.n_solutions), cur_full_paths.n_events)
    cur_events_df['chain'] = cur_events_df['id'] + ':' + cur_events_df['chain_nr'].astype(str)
    cur_events_df['events_per_chrom'] = cur_full_paths.n_events
    cur_events_df['n_paths'] = int(cur_full_paths.n_solutions)
    return cur_events_df


def set_logging_level(logger, level):
    logging_values = {
        'silent': 'WARNING',
        'verbose': 'DEBUG',
    }
    level = logging_values.get(level, level.upper())
    assert level in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    logger.setLevel(level)

 # TODO: move
def get_diffs_from_events_df(cur_id, events_df, supported_chains_only=False):
    cur_events = events_df.query("id == @cur_id")

    # needs to be separate query in case cur_events do not have supported_chain colum
    if supported_chains_only:
        cur_events = events_df.query("supported_chain")
    if cur_events.empty:
        print("no events found")
        return None
    if 'chain_nr' not in cur_events.columns:
        cur_events = cur_events.copy()
        cur_events['chain_nr'] = 0

    all_diffs = []
    for chain_nr, events in cur_events.groupby("chain_nr"):
        diffs = np.stack(
            [
                np.fromiter(diff, int) * (1 if type == "gain" else -1)
                for diff, type in events[["diff", "type"]].values
            ]
        )
        # sort pre before post
        diffs = diffs[np.argsort(events["wgd"].values)[::-1]]
        # revert to get pre-WGD LOHs in the right order
        n_pre = (events["wgd"].values == 'pre').sum()
        diffs[n_pre:] = diffs[n_pre:][::-1]
        all_diffs.append(diffs)
    return all_diffs


def suppress_warnings(warning_type=None):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            with warnings.catch_warnings():
                if warning_type:
                    warnings.simplefilter("ignore", warning_type)
                else:
                    warnings.simplefilter("ignore")
                return func(*args, **kwargs)
        return wrapper
    return decorator


def log_debug(logger, msg):
    if logger.level == logging.DEBUG:
        logger.debug(msg)

def save_fail_reports(failed_reports, results_dir=None, logger=None):
    if results_dir is None:
        from spice import directories, config
        results_dir = directories['results_dir']
    if len(failed_reports) == 0:
        return pd.DataFrame(columns=['id', 'step', 'error', 'status'])
    df_fail = pd.DataFrame(failed_reports)
    for col in ['id', 'step', 'error', 'status']:
        if col not in df_fail.columns:
            df_fail[col] = None
    fail_path = os.path.join(results_dir, config['name'], 'failed_reports.tsv')
    df_fail[['id', 'step', 'error', 'status']].to_csv(fail_path, sep='\t', index=False)
    if logger is not None:
        logger.info(f"A total of {len(df_fail)} tasks failed during execution. "
                    f"Saved failure report with to {fail_path}")
    return df_fail