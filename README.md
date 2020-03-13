# OpenCG

An open-source package for multi-scale coarse-graining (MSCG) modeling in computational chemistry and biology




### Quick Guide

**1. Prerequisites**

Up-to-date Python package management tool is recommended: `conda` or `pip`. 

**2. Download precompiled distributions**

_Python 3.6.x_
```
wget --no-check-certificate https://software.rcc.uchicago.edu/mscg/downloads/opencg-0.0.1-cp36-cp36m-linux_x86_64.whl
```

_Python 3.7.x_
```
wget --no-check-certificate https://software.rcc.uchicago.edu/mscg/downloads/opencg-0.0.1-cp37-cp37m-linux_x86_64.whl
```

_Python 3.8.x_
```
wget --no-check-certificate https://software.rcc.uchicago.edu/mscg/downloads/opencg-0.0.1-cp38-cp38-linux_x86_64.whl
```

**3. Install the distribution**

```
pip install opencg-*.whl --user
```

It is recommended to install the package in the user's base directory, i.e., `~/.local` in the UNIX/Linux environment.

**4. Test the installation**

```
export PATH=${HOME}/.local/bin:${PATH}
cginfo
```
