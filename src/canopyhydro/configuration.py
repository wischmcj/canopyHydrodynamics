from __future__ import annotations

import logging
import os

import toml
import yaml

# Read in environment variables, set defaults if not present

config_file = os.environ.get("CANOPYHYDRO_CONFIG", "./canopyhydro_config.toml")
log_config = os.environ.get("CANOPYHYDRO_LOG_CONFIG", "./logging_config.yml")

with open(log_config) as f:
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

log = logging.getLogger(__name__)

config_file = None
in_flow_grade_lim = None
output_dir = None
input_dir = None
qsm_cols = {}

try:
    with open(config_file) as f:
        config = toml.load(f)
        in_flow_grade_lim = config["model_parameters"]["in_flow_grade_lim"]
        root_dir = config["directories"]["root_dir"]
        output_dir = config["directories"]["output_dir"]
        input_dir = config["directories"]["input_dir"]
        test_input_dir = config["directories"]["test_input_dir"]
        for column in config["qsm"]:
            qsm_cols[column] = config["qsm"][column]
except Exception as e:
    log.error(f"Error loading configuration variables from {config_file}: {e}")
    raise e
