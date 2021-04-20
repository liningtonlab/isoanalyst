# default configuation object for dereplication and import
import json
from pathlib import Path
from typing import Dict, Optional

CONFIG = {
    "colstomatch": ["rettime", "precmz", "precz"],
    "tolerances": {
        "precmz": ["ppm", 10.0],
        "rettime": ["window", 0.03],
        "precz": [None, None],
    },
    "minreps": 3,
}


def load_config(config_file: Path) -> Dict:
    conf = json.load(config_file.open())
    assert all(x in CONFIG.keys() for x in conf.keys())
    return conf


def get_config(
    config_file: Optional[Path],
    minreps: Optional[int] = None,
    mztol: Optional[float] = None,
    rttol: Optional[float] = None,
) -> Dict:
    if config_file is None:
        conf = CONFIG
    else:
        conf = load_config(config_file)
    if minreps is not None:
        conf["minreps"] = minreps
    if mztol is not None:
        conf["tolerances"]["precmz"] = ["ppm", mztol]
    if rttol is not None:
        conf["tolerances"]["rettime"] = ["window", rttol]
    return conf


def get_mz_tol(config):
    return config.get("tolerances", {}).get("precmz", ("ppm", 10.0))[1]