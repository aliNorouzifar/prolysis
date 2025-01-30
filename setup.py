from setuptools import setup, find_packages

setup(
    name="prolysis",
    license="CC BY-NC 4.0",
    version="0.1.1",
    packages=find_packages(),
    install_requires=[
    "automata-lib ==8.4.0",
    "frozendict==2.4.4",
    "networkx==2.8.8",
    "matplotlib==3.10.0",
    "pm4py==2.7.13.1",
    "numpy==1.26.4",
    "pandas==2.2.3",
    "ruptures==1.1.9",
    "scikit_learn==1.6.0",
    "scipy==1.14.1",
    "seaborn==0.13.2",
    "pyemd==1.0.0",
    "redis==5.2.1",
    "numba==0.61.0",
    "pot==0.9.5",
    ],
    author="Ali Norouzifar",
    author_email="ali.norouzifar@pads.rwth-aachen.de",
    description="A Python package for process mining and analysis",
    url="https://github.com/aliNorouzifar/prolysis",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)