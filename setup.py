import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyBindingCurve",
    version="0.2.1",
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
    python_requires='>=3.6',
)

