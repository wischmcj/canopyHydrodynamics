from __future__ import annotations

import logging
import logging.config
import os

import toml
import yaml

cwd = os.getcwd()
print(f"Current working directory: {cwd}")
# Read in environment variables, set defaults if not present
canhydro_location = os.path.dirname(__file__)
print(f"canoPyHydro installed at {canhydro_location}")

env_config_file = os.environ.get(
    "CANOPYHYDRO_CONFIG", f"{canhydro_location}/canopyhydro_config.toml"
)
env_log_config = os.environ.get(
    "CANOPYHYDRO_LOG_CONFIG", f"{canhydro_location}/logging_config.yml"
)

log = logging.getLogger(__name__)

in_flow_grade_lim = None
output_dir = None
input_dir = None
qsm_cols = {}


def load_config(config_file, log_config):
    try:
        with open(log_config) as f:
            config = yaml.safe_load(f.read())
            logging.config.dictConfig(config)
    except FileNotFoundError:
        log.error("Error loading default log config")
        log_config = f"{canhydro_location}/logging_config.yml"

    log = logging.getLogger(__name__)

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
    except:
        log.error("Error loading default log config")
        config_file = f"{canhydro_location}/logging_config.yml"

    return config_file, log_config


load_status = load_config(env_config_file, env_log_config)

if load_status != (env_config_file, env_log_config):
    load_config(env_config_file, env_log_config)
