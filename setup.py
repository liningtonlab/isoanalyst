from setuptools import setup
from isoanalyst import VERSION

setup(
    name="isoanalyst",
    version=VERSION,
    packages=["isoanalyst", "isoanalyst.convert"],
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
