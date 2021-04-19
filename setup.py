from setuptools import setup

setup(
    name="isoanalyst",
    version="0.1.0",
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
