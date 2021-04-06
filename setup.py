import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pybindingcurve",
    version="1.0.5",
    author="Steven Shave",
    author_email="steve.nshave@gmail.com",
    description="Protein ligand binding simulation in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/stevenshave/pybindingcurve",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        # Excluding numpy 1.19.4 as a bug in the windows 10 update 2004
        # causes a problem - see
        # https://tinyurl.com/y3dm3h86, and
        # https://github.com/numpy/numpy/issues/16744.
        "numpy>=1.18.0,!=1.19.4",
        "matplotlib>=3.2.1",
        "lmfit>=1.0.0",
        "mpmath>=1.1.0",
        "autograd>=1.3",
    ],
    python_requires=">=3.6",
)
