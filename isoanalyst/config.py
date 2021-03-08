# default configuation object for dereplication and import

CONFIG = {
    "ColsToMatch": ["RetTime", "PrecMz", "PrecZ"],
    "Tolerances": {
        "PrecMz": ["ppm", 10],
        "RetTime": ["window", 0.03],
        "PrecZ": [None, None],
    },
    "MinReps": 3,
}
