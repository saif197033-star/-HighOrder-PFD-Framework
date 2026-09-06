from setuptools import setup, find_packages

setup(
    name="saif_raad_pfd",
    version="1.0.0",
    author="Saif R. Lazim",
    author_email="saifraad.math@gmail.com",
    description="A Unified Non-Recursive Framework for High-Order Partial Fraction Decomposition",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/saifraad/saif_raad_pfd",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "sympy>=1.8",
    ],
)