from __future__ import annotations

import time
import logging
import toml
import calendar
from pathlib import Path


with open('canhydro.toml', 'r') as f:
    config = toml.load(f)

DIR = config['directories']['root']

current_GMT = time.gmtime()
time_stamp = str(calendar.timegm(current_GMT))

log_dir = DIR
log_dir = Path(''.join([log_dir,'log',str(time_stamp)]))


logging.basicConfig(filename=log_dir  , filemode='w', level=logging.DEBUG, encoding='utf-8')
LOGGER = logging.getLogger("my-logger")

