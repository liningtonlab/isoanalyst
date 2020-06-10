from setuptools import setup

setup(
name="isotracer",
version="0.0.2",
packages=["isotracer"],
    entry_points={
    'console_scripts': [
        'isotracer = isotracer.cli:main',
    ]
    },
python_requires='>=3.6'
)

