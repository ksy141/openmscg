# OpenCG

An open-source package for multi-scale coarse-graining (MSCG) modeling in computational chemistry and biology




### Quick Guide

**1. Prerequisites**

Up-to-date Python package management tool is recommended: `conda` or `pip`. 

**2. Clone the pre-compiled distribution from GitLab**

```
git clone git@software.rcc.uchicago.edu:MSCG/opencg.git
```

**3. Install the distribution**

```
cd opencg
python setup.py install --user
```

It is recommended to install the package in the user's base directory, i.e., `~/.local` in the UNIX/Linux environment.

**4. Test the installation**

```
export PATH=${HOME}/.local/bin:${PATH}
cgib --help
```
