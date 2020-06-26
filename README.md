# OpenMSCG

An open-source package for multi-scale coarse-graining (MSCG) modeling in computational chemistry and biology




### Quick Guide

**1. Prerequisites**

Up-to-date Python package management tool is recommended: `conda` or `pip`. 

**2. Download precompiled distributions**

_For Python 3.6.x users:_
```
wget --no-check-certificate https://software.rcc.uchicago.edu/download/OpenMSCG-0.1.0-cp36-cp36m-linux_x86_64.whl
```

_For Python 3.7.x users:_
```
wget --no-check-certificate https://software.rcc.uchicago.edu/download/OpenMSCG-0.1.0-cp37-cp37m-linux_x86_64.whl
```

_For Python 3.8.x users:_
```
wget --no-check-certificate https://software.rcc.uchicago.edu/download/OpenMSCG-0.1.0-cp38-cp38-linux_x86_64.whl
```

**3. Install the distribution**

```
pip install OpenMSCG-*.whl --user
```

It is recommended to install the package in the user's base directory, i.e., `~/.local` in the UNIX/Linux environment.

**4. Test the installation**

```
export PATH=${HOME}/.local/bin:${PATH}
cginfo
```


