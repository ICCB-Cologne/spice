import os
from pathlib import Path
import yaml
import collections.abc
from typing import Optional
import sys

# Use importlib.resources for accessing package data (works with installed packages)
if sys.version_info >= (3, 9):
    from importlib.resources import files
else:
    try:
        from importlib_resources import files
    except ImportError:
        # Fallback for older Python versions without importlib_resources
        files = None

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

def _read_yaml_content(content):
    """Read YAML from string content."""
    return yaml.safe_load(content)

def _read_yaml(path):
    """Read YAML from a filesystem path."""
    with open(path, 'rt') as f:
        return yaml.safe_load(f.read())

def _get_default_config_content():
    """Get default_config.yaml content via package resources."""
    if files is None:
        raise FileNotFoundError("importlib.resources unavailable for default_config.yaml")
    try:
        return files('spice').joinpath('objects', 'default_config.yaml').read_text()
    except (TypeError, ImportError, AttributeError, FileNotFoundError) as exc:
        raise FileNotFoundError("Could not find default_config.yaml in spice/objects/") from exc

def load_config(config_path: Optional[str] = None, assert_exists: bool = True):
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
    default_config_content = _get_default_config_content()
    default_config = _read_yaml_content(default_config_content)

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

    if config_path is not None:
        # Require base_dir and derive fallbacks from it
        base_dir = directories.get('base_dir')
        if not base_dir:
            raise ValueError("Config file must specify 'directories.base_dir'.")
        base_dir_path = Path(base_dir).expanduser().resolve()
        directories['base_dir'] = str(base_dir_path)
        fallback_dirs = {
            'data_dir': base_dir_path / 'data',
            'results_dir': base_dir_path / 'results',
            'log_dir': base_dir_path / 'logs',
            'plot_dir': base_dir_path / 'plots',
            'tmp_dir': base_dir_path / 'tmp',
        }
        for k, v in fallback_dirs.items():
            if not directories.get(k):
                directories[k] = str(v)

        # Persist back into config
        config['directories'] = directories

        # Normalize all directory paths: if not absolute, resolve under base_dir
        for dkey, dval in list(directories.items()):
            if isinstance(dval, str) and dval and not os.path.isabs(dval):
                directories[dkey] = str((base_dir_path / dval).resolve())
        config['directories'] = directories

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
            cur_abs = os.path.abspath(os.environ['SPICE_CONFIG'])
            if cur_abs != new_abs:
                import warnings
                warnings.warn(
                    f"SPICE_CONFIG already set to '{cur_abs}'. Requested different path '{new_abs}' will be ignored.")
