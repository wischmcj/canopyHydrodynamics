from __future__ import annotations

import calendar
import logging
import sys
import time
from pathlib import Path

import toml
from memory_profiler import LogFile

with open("user_def_config.toml") as f:
    config = toml.load(f)

# Load Dirs
DIR = config["directories"]["root"]

input_dir = DIR
input_dir = Path("".join([input_dir, "data", "\\input"]))

output_dir = DIR
output_dir = Path("".join([output_dir, "data", "\\output"]))

test_input_dir = DIR
test_input_dir = Path("".join([test_input_dir, "data", "\\test"]))

# Current datetime
current_GMT = time.gmtime()
time_stamp = str(calendar.timegm(current_GMT))

# QSM column order
qsm_cols = {}
for column in config["qsm"]:
    qsm_cols[column] = config["qsm"][column]

config_vars = {}
for column in config["config_vars"]:
    config_vars[column] = config["config_vars"][column]

log_dir = DIR
log_dir = Path("".join([log_dir, r"log\log_", str(time_stamp)]))

sys.stdout = LogFile(str(log_dir))

logging.basicConfig(
    filename=log_dir, filemode="w", level=logging.INFO, encoding="utf-8"
)
LOGGER = logging.getLogger("my-logger")
log = logging.getLogger("my-logger")
