Introduction
============

The multiscale coarse-graining (MSCG) methodology provides a systematic,
bottom-up way to calculate effective CG interactions based on rigorous
statistical mechanics. It seeks to approximate the many-body potential of
mean force by variationally minimizing the difference between CG forces and
the mapped fine-grained forces (a.k.a, "force matching"). This
method is related to liquid state theory and the Yvon-Born-Green equation.

This project is funded by NSF/SSE (No.5-27425), and focused on the sustainable
development of a cutting-edge CG/UCG molecular simulation platform with broad
applicability to the fields of chemistry, biology, and materials science
in accordance with the National Strategic Supercomputing Initiative (e.g.
Principle 1, Objective 2), and the Vision and Strategy for Software for
Science, Engineering, and Education (NSF 12-113, Strategies 1-4).

The OpenMSCG package provides an implementation of the algorithms described in the following
references. Reading these papers before using the code will facilitate understanding of the
terminology throughout this manual and the comments throughout the source code.

.. seealso::
    L. Lu, S. Izvekov, A. Das, H. C. Andersen, and G. A. Voth, "Efficient, Regularized, and
    Scalable Algorithms for Multiscale Coarse-Graining", Journal of Chemical Theory and
    Computation, 6(3), 954-965 (2010). doi:10.1021/ct900643r

**This package also provides support for the methods described in these references:**

* "Multiscale coarse graining of liquid-state systems", S. Izvekov, and G. A. Voth, The Journal of Chemical Physics, 123, 134105 (2005). doi:10.1063/1.2038787

* "The multiscale coarse-graining method. II. Numerical implementation for coarse-grained molecular models," W. G. Noid, P. Liu, Y. Wang, J.-W. Chu, G. S. Ayton, S. Izvekov, H. C. Andersen, and G. A. Voth, The Journal of Chemical Physics, 128, 244115 (2008). doi:10.1063/1.2938857

* "The theory of ultra-coarse-graining. 2. Numerical implementation" A. Davtyan, J. F. Dama, A. V. Sinitskiy, and G. A. Voth, The Journal of Chemical Theory and Computation, 10, 5265–5275 (2014). doi:10.1021/ct500834t

* "The theory of ultra-coarse-graining. 3. Coarse-grained sites with rapid local equilibrium of internal states", J. F. Dama, J. Jin, and G. A. Voth, The Journal of Chemical Theory and Computation, 13, 1010-1022 (2017). doi: 10.1021/acs.jctc.6b01081

* "Coarse-graining errors and numerical optimization using a relative entropy framework," A. Chaimovich and M. Shell, The Journal of Chemical Physics, 134, 094112 (2011). doi:10.1063/1.3557038
