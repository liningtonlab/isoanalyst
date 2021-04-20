from setuptools import setup
from isoanalyst import VERSION

setup(
    name="isoanalyst",
    version=VERSION,
    packages=["isoanalyst"],
    install_requires=[
        "click",
    ],
    entry_points={
        "console_scripts": [
            "isoanalyst = isoanalyst.cli:cli",
        ]
    },
    python_requires=">=3.6",
)
