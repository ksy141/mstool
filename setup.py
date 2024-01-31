"""
Installation script for mstool

This assumes that the python (>=3.7) that is used to execute this script
is the conda / miniconda environment into which you want to install mstool.

USAGE:

    python -m build

or 

    pip install .

"""

#import sys
#import os
#import subprocess
#import shutil
#
## We can only run this script from the sire directory
#curdir = os.path.abspath(".")
#
## Find the conda base
#conda_base = os.path.abspath(os.path.dirname(sys.executable))
#if os.path.basename(conda_base) == "bin":
#    conda_base = os.path.dirname(conda_base)
#
## MacOS and Linux
#conda_bin  = os.path.join(conda_base, "bin")
#python_bin = os.path.join(conda_bin, "python")
#conda = os.path.join(conda_bin, "conda")

# Read version
def _read_version(ifile):
    with open(ifile) as fp:
        for line in fp:
            line = line.strip()
            if not line: continue
            sl = [s.strip() for s in line.split('=')]
            if sl[0] == '__version__':
                return sl[1].strip('\"')
    
# Find the installed dependencies
def _get_installed():
    installed = {}
    p = subprocess.run(["conda", "list"], stdout=subprocess.PIPE)
    lines = p.stdout.decode().split('\n')
    for line in lines:
        line = line.strip()
        line = line.split('#')[0]
        if not line: continue
        sl = line.split()
        installed[sl[0]] = sl[1]
    return installed

# Get Requirements
def _get_required():
    import re
    required = {}
    with open('requirements.txt') as fp:
        for line in fp:
            line = line.strip()
            if not line: continue
            words = [word.strip() for word in re.split("[<>]*=", line)]
            if len(words) == 1:
                required[words[0]] = '0.0'
            else:
                required[words[0]] = words[1]
    return required

# conda install
def _conda_install(dependency):
    cmd = f'{conda} install -c conda-forge {dependency}'
    print(cmd)
    subprocess.run(cmd.split())
            
def _build():
    os.chdir("src/mstool/lib")
    subprocess.run(["python", "setup.py", "build_ext", "--inplace"])
    os.chdir(curdir)

def _compare_version(str1, str2):
    str1s = str1.split('.')
    str2s = str2.split('.')

    if len(str1s) >= len(str2s):
        t = len(str2s)
    else:
        t = len(str1s)

    for i in range(t):
        if int(str1s[i]) > int(str2s[i]):
            return 1
        if int(str1s[i]) < int(str2s[i]):
            return 2


### INSTALL
#installed = _get_installed()
#required  = _get_required()

#for pkg, req_ver in required.items():
#    if pkg not in installed.keys():
#        _conda_install(pkg)
#    else:
#        ins_ver = installed[pkg]
#        com = _compare_version(req_ver, installed[pkg])
#
#        if com == 1:
#            print(f"{pkg} is already installed but needs to be updated")
#            _conda_install(pkg)
#
#        if com == 2:
#            print(f"{pkg} {ins_ver} is already installed and fine") 
 

### Build
#_build()

# why do I have to specify packages in setup and in pyproject.toml
# [tool.setuptools.packages.find]
# where = ["src/mstool"]
#    version=_read_version('src/mstool/version.py'),

from setuptools   import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np


ext_modules=[
    Extension("mstool.lib.distancelib",
             ["src/mstool/lib/distancelib.pyx"]),
    Extension("mstool.lib.qcprot",
             ["src/mstool/lib/qcprot.pyx"]) 
]


setup(
    version=_read_version('src/mstool/version.py'),
    package_dir={"": "src"},
    include_package_data = False,
    package_data={
        "mstool.lib": ["*.pyx"],
        "mstool.examples.Backmapping": ["*/cg*pdb", "*/protein_AA.pdb", "*/*.xml", "*/*.dat"],
        "mstool.FF": ["*/*.xml", "*/*.itp", "*/*.pdb"],
        "mstool.mapping": ["*.dat"],
        "mstool.mapping.structures": ["*.pdb"],
    },
    exclude_package_data={
        "mstool.backup": ["*.py"]
    },
    ext_modules = ext_modules,
    include_dirs = [np.get_include()],
)

