from __future__ import annotations

import calendar
import logging
import time
from pathlib import Path

import toml

with open("user_def_config.toml") as f:
    config = toml.load(f)

#Load Dirs
DIR = config["directories"]["root"]

input_dir=DIR
input_dir = Path("".join([input_dir,"data", "input"]))

output_dir = DIR
output_dir = Path("".join([output_dir,"data", "output"]))

#Current datetime
current_GMT = time.gmtime()
time_stamp = str(calendar.timegm(current_GMT))

#QSM column order
qsm_cols = {}
for column in config["qsm"]:
    qsm_cols[column] = config["qsm"][column]
print(qsm_cols)

log_dir = DIR
log_dir = Path("".join([log_dir, "log", str(time_stamp)]))


logging.basicConfig(
    filename=log_dir, filemode="w", level=logging.DEBUG, encoding="utf-8"
)
LOGGER = logging.getLogger("my-logger")
