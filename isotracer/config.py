# default configuation object for dereplication
# could consider having all config options here

CONFIG = {
    "ColsToMatch": ["RetTime", "PrecMz", "PrecZ"],
    "Tolerances": {
        "PrecMz": ["ppm", 10],
        "RetTime": ["window", 0.03],
        "PrecZ": [None, None],
    },
    "MinReps": 3,
}
