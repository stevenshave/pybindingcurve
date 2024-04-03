import io
import os
import re

from setuptools import find_packages
from setuptools import setup


def read(filename):
    filename = os.path.join(os.path.dirname(__file__), filename)
    text_type = type(u"")
    with io.open(filename, mode="r", encoding="utf-8") as fd:
        return re.sub(text_type(r":[a-z]+:`~?(.*?)`"), text_type(r"``\1``"), fd.read())


from distutils.core import setup

setup(
    name="PyBindingCurve",
    version="1.1.0",
    packages=find_packages(exclude=("tests",)),
    license="The MIT License (MIT)",
    long_description=open("README.md").read(),
)
