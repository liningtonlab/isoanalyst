from setuptools import setup

setup(
name="isoanalyst",
version="0.0.2",
packages=["isoanalyst"],
    entry_points={
    'console_scripts': [
        'isoanalyst = isoanalyst.cli:main',
    ]
    },
python_requires='>=3.6'
)

