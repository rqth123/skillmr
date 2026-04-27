# setup.py
from setuptools import setup, find_packages

setup(
    name="gwas-agent",
    version="1.0.0",
    description="GWAS Analysis Agent for OpenClaw",
    author="OpenClaw Team",
    packages=find_packages(),
    install_requires=[
        "openclaw>=0.1.0",
        "requests>=2.28.0",
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "matplotlib>=3.6.0",
        "scipy>=1.9.0",
    ],
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
