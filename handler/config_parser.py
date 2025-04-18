import yaml # type: ignore
import os
from typing import Any, Dict, List, Optional, Union

# Default config file path
DEFAULT_CONFIG_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config.yaml")

# Default configuration template
DEFAULT_CONFIG = {
    "data_root_path": "/big8_disk/data",
    "run_root_path": os.path.join("/bip8_disk", os.getenv("USER")),
    "ans_rank": {
        "NYGC": 2,
        "SEQC2": 2,
        "orthogonal-tools-benchmark": 1
    },
    "ans_vcf_filter": [
        "superSet"
    ],
    "sequencing_platform_rank": {
        "ONT_5khz_simplex_5mCG_5hmCG": -1,
        "ONT": 1,
        "ONT_Dorado": -1,
        "ONT_PAO": 1
    },
    "run_samples": [
        "H2009",
        "HCC1954",
        "HCC1937",
        "COLO829",
        "H1437",
        "HCC1395"
    ],
    "run_puritys": [
        1,
        0.8,
        0.6,
        0.4,
        0.2
    ],
    "run_calling_softwares": [
        "DeepSomatic_TO_v1_8_0",
        "ClairS_TO_v0_3_0",
        "ClairS_TO_v0_3_0_pileup",
        "ClairS_TO_v0_3_0_pileup_nonsomatic",
        "Longphase_TO_v0_0_1"
    ],
    "run_type": [
        "snv"
    ]
}

def generate_default_config(output_path: Optional[str] = None, overwrite: bool = False) -> str:
    """
    Generate a default configuration file.
    
    Args:
        output_path: Path where to save the config file. If None, uses the default path.
        overwrite: Whether to overwrite the file if it already exists.
        
    Returns:
        Path to the generated config file.
        
    Raises:
        FileExistsError: If the file already exists and overwrite is False.
        PermissionError: If the file can't be written to.
    """
    if output_path is None:
        output_path = DEFAULT_CONFIG_PATH
    
    # Check if file already exists
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(f"Config file already exists at: {output_path}. Use overwrite=True to replace it.")
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Write the default config
    try:
        with open(output_path, "w") as f:
            yaml.dump(DEFAULT_CONFIG, f, default_flow_style=False)
        return output_path
    except PermissionError:
        raise PermissionError(f"Cannot write to config file: {output_path}")

def get_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load and return the configuration from a YAML file.
    
    Args:
        config_path: Path to the config file. If None, uses the default path.
        
    Returns:
        Dict containing the configuration.
        
    Raises:
        FileNotFoundError: If the config file doesn't exist.
        yaml.YAMLError: If the YAML file is invalid.
    """
    if config_path is None:
        config_path = DEFAULT_CONFIG_PATH
    
    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
        return config
    except FileNotFoundError:
        raise FileNotFoundError(f"Config file not found at: {config_path}")
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"Error parsing YAML file: {e}")

def get_config_value(key_path: Union[str, List[str]], default: Any = None, config_path: Optional[str] = None) -> Any:
    """
    Get a value from the config using a dot-notation path or list of keys.
    
    Args:
        key_path: String with dot notation (e.g., "ans_rank.NYGC") or list of keys ["ans_rank", "NYGC"]
        default: Default value to return if the key doesn't exist
        config_path: Optional path to the config file
        
    Returns:
        The value at the specified path or the default value if not found
    """
    config = get_config(config_path)
    
    # Convert string path to list of keys
    if isinstance(key_path, str):
        keys = key_path.split('.')
    else:
        keys = key_path
    
    # Navigate through the nested dictionary
    current = config
    try:
        for key in keys:
            current = current[key]
        return current
    except (KeyError, TypeError):
        return default

def get_config_column(column_name: str, config_path: Optional[str] = None) -> Any:
    """
    Get a top-level column from the config.
    
    Args:
        column_name: Name of the column to retrieve
        config_path: Optional path to the config file
        
    Returns:
        The value of the specified column
        
    Raises:
        KeyError: If the column doesn't exist
    """
    config = get_config(config_path)
    try:
        return config[column_name]
    except KeyError:
        raise KeyError(f"Column '{column_name}' not found in config")

def get_my_address(config_path: Optional[str] = None) -> str:
    """
    Get the user's computer address in the format user@ip.
    
    Args:
        config_path: Optional path to the config file
        
    Returns:
        String in the format "user@ip"
    """
    config = get_config(config_path)
    try:
        return f"{config['my_computer']['user']}@{config['my_computer']['ip']}"
    except KeyError:
        raise KeyError("Missing 'my_computer' configuration with 'user' and 'ip' fields")

def get_boos_address(config_path: Optional[str] = None) -> str:
    """
    Get the boos address in the format user@ip.
    
    Args:
        config_path: Optional path to the config file
        
    Returns:
        String in the format "user@ip"
    """
    config = get_config(config_path)
    try:
        return f"{config['boos']['user']}@{config['boos']['ip']}"
    except KeyError:
        raise KeyError("Missing 'boos' configuration with 'user' and 'ip' fields")

def get_boos_default_path(config_path: Optional[str] = None) -> str:
    """
    Get the boos default path.
    
    Args:
        config_path: Optional path to the config file
        
    Returns:
        The default path string
    """
    config = get_config(config_path)
    try:
        return config['boos']['default_path']
    except KeyError:
        raise KeyError("Missing 'boos.default_path' configuration")

def update_config(updates: Dict[str, Any], config_path: Optional[str] = None) -> None:
    """
    Update the configuration file with new values.
    
    Args:
        updates: Dictionary of updates to apply to the config
        config_path: Optional path to the config file
        
    Raises:
        FileNotFoundError: If the config file doesn't exist
        PermissionError: If the file can't be written to
    """
    if config_path is None:
        config_path = DEFAULT_CONFIG_PATH
    
    # Read current config
    config = get_config(config_path)
    
    # Apply updates
    def deep_update(d, u):
        for k, v in u.items():
            if isinstance(v, dict) and k in d and isinstance(d[k], dict):
                deep_update(d[k], v)
            else:
                d[k] = v
    
    deep_update(config, updates)
    
    # Write back to file
    try:
        with open(config_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)
    except PermissionError:
        raise PermissionError(f"Cannot write to config file: {config_path}")
