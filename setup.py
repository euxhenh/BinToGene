import setuptools
from setuptools import Command
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = [
    'numpy',
    'pandas'
]

class CleanCommand(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

cmdclass = {'clean': CleanCommand}

setuptools.setup(
    name="BinToGene",
    package_dir={'BinToGene': 'src'},
    packages=['BinToGene'],
    version="1.0",
    author="Euxhen Hasanaj",
    author_email="ehasanaj@cs.cmu.edu",
    description="A bin to gene conversion package.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ferrocactus/BinToGene",
    install_requires=install_requires,
    cmdclass=cmdclass,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
