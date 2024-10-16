from __future__ import annotations

import logging
import logging.config

import yaml

with open("logging_config.yml") as f:
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

script_logger = logging.getLogger("script")
