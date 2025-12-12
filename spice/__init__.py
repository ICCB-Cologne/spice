import os
from pathlib import Path
import yaml
import collections.abc

def update_nested_dict(orig_dict, updated_values):
    '''
    See https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    '''
    orig_dict = orig_dict.copy()
    for k, v in updated_values.items():
        if isinstance(v, collections.abc.Mapping):
            orig_dict[k] = update_nested_dict(orig_dict.get(k, {}), v)
        else:
            orig_dict[k] = v
    return orig_dict

# Module-level globals that can be updated by CLI before submodules import.
default_config = None
config = None
directories = None

def _read_yaml(path):
    with open(path, 'rt') as f:
        return yaml.safe_load(f.read())

def load_config(config_path: str | None = None, assert_exists: bool = True):
    """
    Load SPICE configuration by merging default_config.yaml with an optional user config.

    If config_path is None, attempts to read a workspace-level 'config.yaml' next to
    default_config.yaml. This function updates module-level globals: `config` and `directories`.
    """
    global default_config, config, directories

    if config_path is not None:
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found at provided path: '{config_path}'.")
        set_config(config_path)

    # Always load the default first
    default_path = os.path.join(Path(__file__).parent.parent, 'default_config.yaml')
    default_config = _read_yaml(default_path)

    # Determine user config source
    user_cfg = {}
    # Environment variable override
    env_cfg_path = os.environ.get('SPICE_CONFIG')
    cfg_candidate = config_path or env_cfg_path
    if cfg_candidate is not None and os.path.exists(cfg_candidate):
        user_cfg = _read_yaml(cfg_candidate) or {}

    # Merge user config onto defaults
    config = update_nested_dict(default_config, user_cfg)
    # Track whether we're still using the default-only configuration
    if isinstance(config, dict):
        config.setdefault('meta', {})
        config['meta']['is_default'] = (user_cfg == {} and env_cfg_path is None and config_path is None)
        if config_path or env_cfg_path:
            # Record where overrides came from
            config['meta']['source_path'] = os.path.abspath(config_path or env_cfg_path)
    directories = config.get('directories', {})

    # Fallbacks for core directories if not provided
    package_root = Path(__file__).parent.parent
    directories['base_dir'] = str(package_root)
    fallback_dirs = {
        'data_dir': package_root / 'data',
        'results_dir': package_root / 'results',
        'log_dir': package_root / 'logs',
        'plot_dir': package_root / 'plots',
        'tmp_dir': package_root / 'tmp',
    }
    for k, v in fallback_dirs.items():
        if not directories.get(k):
            directories[k] = str(v)

    # Persist back into config
    config['directories'] = directories

    # Normalize all directory paths: if not absolute, resolve under package root
    for dkey, dval in list(directories.items()):
        if isinstance(dval, str) and dval and not os.path.isabs(dval):
            directories[dkey] = str((package_root / dval).resolve())
    config['directories'] = directories

    if 'input_files' in config and 'knn_train' in config['input_files']:
        knn_path = config['input_files']['knn_train']
        if knn_path and not os.path.isabs(knn_path):
            abs_knn_path = str((package_root / knn_path).resolve())
            config['input_files']['knn_train'] = abs_knn_path
    return config

# Perform an initial load so library usage without CLI still works as before
load_config()

def set_config(config_path=None):
    """Set the SPICE_CONFIG environment variable for worker processes.

    Synchronize environment for worker processes with safeguards
    Only set if a valid explicit config_path was provided and "SPICE_CONFIG" wasn't already set.
    """

    if config_path and os.path.exists(config_path):
        new_abs = os.path.abspath(config_path)
        if not os.environ.get('SPICE_CONFIG'):
            os.environ['SPICE_CONFIG'] = new_abs
        else:
            cur_abs = os.path.abspath(cur_env)
            if cur_abs != new_abs:
                import warnings
                warnings.warn(
                    f"SPICE_CONFIG already set to '{cur_abs}'. Requested different path '{new_abs}' will be ignored.")
