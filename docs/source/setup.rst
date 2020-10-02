Compilation
=================

This document details instructions to download, compile, and install the OpenMSCG package from the source code.
If you're not interested in optimizing or developing the package, it is recommended that you use the 
easy-install method with our precompiled binary distributions in wheel format.

Setup SSH Key on GitLab
-----------------------

Login into the GitLab portal (`Link <https://software.rcc.uchicago.edu/git/>`_)
from the RCC Software website. If you don't have a GitLab account, you can sign up 
at this `page <https://software.rcc.uchicago.edu/mscg/register.php>`_ for free.

Setup your SSH public key to connect to GitLab (`instructions <https://help.github.
com/en/articles/adding-a-new-ssh-key-to-your-github-account>`_ ). If you are unfamiliar
with SSH key pairs or how to create them, please read the instructions
`here <https://www.digitalocean.com/docs/droplets/how-to/add-ssh-keys/create-
with-openssh/>`_.

Check that the SSH key was paired successfully::

    ssh -T git@software.rcc.uchicago.edu
    Welcome to GitLab, @yuxing!

If the SSH key is working, a welcome message from GitLab will be shown.


Requirements
------------

To compile the OpenMSCG package you need to have the following software packages 
in your computer system:

1. Python (version>=3.6)
2. C++ compiler (support C++11 standard)
3. GSL library (version>=2.1)
4. LAPACK library

**IMPORTANT**: Python 2 is not supported.

**IMPORTANT**: If you're working on HPC systems, these packages have likely been 
installed by the system administrators already. Please read the system manual 
or contact your support staff for instructions on how to use them.


1. Python3
^^^^^^^^^^

There are three installation options for the Python libraries needed by the
package. We recommend using Anaconda3.

**Option 1: Anaconda3**

``Anaconda3`` can be downloaded and installed from the website::

    https://www.anaconda.com/distribution/#download-section

Please choose the distribution with ``Python 3.7 version``

**Option 2: Miniconda3**

As indicated by its name, ``Miniconda`` is a mini-version of ``Anaconda``, 
which only comes with a small number of packages to minimize its disk usage. 
``Miniconda`` also provides very friendly functions for package 
management::

    https://docs.conda.io/en/latest/miniconda.html

Please choose the distribution that fits your OS and with ``Python 3.7``.

**Option 3: From official Python3 distribution:**

The software can be downloaded from the website::

    https://www.python.org/downloads/

But the installation requires more effort.

.. note:: If you are using Anaconda, you do not need to worry about installing anything else and simply use the quick start guide. If not, continue.


2. C++ Compiler
^^^^^^^^^^^^^^^

**Option 1: GNU Compiler (g++)**

On HPC Linux systems, GCC is always provided. Please check with the system 
manual or contact assistance to initialize the GCC environment. If you're 
using Linux or MacOS on a local machine, you 
can follow the instructions below to install it:

* `Ubuntu <https://linuxconfig.org/how-to-install-gcc-the-c-compiler-on-ubuntu-18-04-bionic-beaver-linux>`_
* `Redhat/CentOS <https://www.cyberciti.biz/faq/centos-rhel-7-redhat-linux-install-gcc-compiler-development-tools/>`_
* `MacOS <https://www.cyberciti.biz/faq/howto-apple-mac-os-x-install-gcc-compiler/>`_

**Note:** Version 4.3 or above is required to support C++11.

**Option 2: Intel Compiler (icpc)**

Intel compiler usually gives better performance on an Intel processor, 
however, it's not free. You need to buy licenses to install and use it. 
On HPC Linux systems using Intel processors the Intel compiler is usually 
provided. Please check with the system manual or contact assistance 
to initialize the ICC/ICPC environment.

3. GSL Library
^^^^^^^^^^^^^^

The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers. 
It is a free open source library under the GNU General Public License. By default, 
you don't need to compile the library by yourself. The GNU website provides a list of 
compiled binaries to download. You need to choose the correct platform and 
version. A general guide can be found at::

    https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/

Some simple commands on common OS's to automatically install GSL:

* **Ubuntu 18.04 or above**: ``sudo apt-get install libgsl-dev``

* **Redhat/CentOS 7 or above**: ``sudo yum install gsl-devel``

* **MacOS X or above**: ``brew install gsl`` (`GSL website <http://macappstore.org/gsl/>`_)


4. LAPACK Library
^^^^^^^^^^^^^^^^^

Quoted from the `LAPACK <http://www.netlib.org/lapack/>`_ official website: "LAPACK
is written in Fortran 90 and provides routines for solving systems of 
simultaneous linear equations, least-squares solutions of linear systems of 
equations, eigenvalue problems, and singular value problems."

On most Linux systems, LAPACK is already installed as a built-in library.
Therefore, there's often no need to manually set it up. In case it is not 
found, you can follow the instructions from the following website to install 
it::

    https://coral.ise.lehigh.edu/jild13/2016/07/27/install-lapack-and-blas-on-linux-based-systems/

**Setup with Intel MKL:**

If you're using Intel compilers, using Intel's Math Kernel Library (MKL) to provide
LAPACK interfaces can give much better performance (note: this is only recommended 
if you are using Intel processors). MKL is free to be obtained with instructions on 
the following website::

    https://software.intel.com/en-us/get-started-with-mkl-for-linux

Again, on HPC systems, MKL should be provided, and is usually found with the Intel 
compilers.

**Setup on MacOS:**

Instead of the standard LAPACK library, Apple provides its own framework ``Accelarate`` 
to use, which can be used by passing the flag ``-framework Accelerate`` to the linker
command.


Download Source
---------------

In the folder that you want to download the package, run "git" to clone the remote
repository::

    git clone git@software.rcc.uchicago.edu:MSCG/OpenMSCG.git

Enter the repository directory::

    cd OpenMSCG


Configure
---------

To compile the C++ code for the package, you need to prepare a configuration file named
as ``build.cfg`` in the root folder to specify required options. An example file is prepared
in the ``arch`` folder with the following content ::

    [build_options]
    cc         = icpc
    cxx        = icpc
    compile    = -O2 -Wno-sign-compare
    link       = -static -static-intel -static-libgcc -static-libstdc++ -wd10237
    gsl_lib    = /software/gsl-2.2.1-el7-x86_64+intel-16.0/lib/libgsl.a
                 /software/gsl-2.2.1-el7-x86_64+intel-16.0/lib/libgslcblas.a
    lapack_lib = -mkl=sequential -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -liomp5


These options may need to be customized to match your system envrionment:

  * ``cc``: command name for C++ compilation
  * ``cxx``: command name for C++ compilation, the same as ``cc``
  * ``compile``: options for compilation, for examples, extra header files for GSL and Linpack.
  * ``link``: options for building (linking) target binaries
  * ``gsl_lib``: options specifying the GSL library for linking
  * ``lapack_lib``: options specifying the Linpack libraries for linking



Build and Install
-----------------

The package is prepared using Python 
`SetupTools <https://setuptools.readthedocs.io/en/latest/setuptools.html>`_ . 
Run the following command ::

    python setup.py build_ext --inplace
    python setup.py install

If you wish to install it in your local home directory, you can include the ``--user`` flag.



















