Introduction
============

The multiscale coarse-graining (MSCG) methodology provides a systematic,
bottom-up way to calculate effective CG interactions based on rigorous
statistical mechanics. It seeks to approximate the many-body potential of
mean force by variationally minimizing the difference between CG forces at
the mapped fine-grained reference forces (a.k.a, "force matching"). This
method is related to liquid state theory and the Yvon-Born-Green equation.

This project is funded by NSF/SSE (No.5-27425), and focused on the sustainable
development of a cutting-edge CG/UCG molecular simulation platform with broad
applicability to the fields of chemistry, biology and the materials sciences
in accordance with the National Strategic Supercomputing Initiative (e.g.
Principle 1, Objective 2), and the Vision and Strategy for Software for
Science, Engineering, and Education (NSF 12-113, Strategies 1-4).

The OpenMSCG package provides an implementation of the algorithms described in this
reference. Please read this paper before using the code to understand the terminology of
this manual and the comments in the source code.

.. seealso::
    L. Lu, S. Izvekov, A. Das, H. C. Andersen, and G. A. Voth, "Efficient, Regularized, and
    Scalable Algorithms for Multiscale Coarse-Graining", Journal of Chemical Theory and
    Computation, 6(3), 954-965 (2010). doi:10.1021/ct900643r

**This package also provides support for the methods described in these references:**

* "Multiscale coarse graining of liquid-state systems", S. Izvekov, and G. A. Voth, The Journal of Chemical Physics, 123(13), 134105 (2005). doi:10.1063/1.2038787

* "A Bayesian statistics approach to multiscale coarse graining", P. Liu, Q. Shi, H. D. Daume III, and G. A. Voth, Journal of Chemical Physics, 129, 214114 (2008). doi:10.1063/1.3033218

* "Multiscale coarse graining of liquid-state systems", S. Izvekov, and G. A. Voth, The Journal of Chemical Physics, 123(13), 134105 (2005). doi:10.1063/1.2038787

* "A Bayesian statistics approach to multiscale coarse graining", P. Liu, Q. Shi, H. D. Daume III, and G. A. Voth, Journal of Chemical Physics, 129, 214114 (2008). doi:10.1063/1.3033218

* "The multiscale coarse-graining method. VI. Implementation of three-body coarse-grained potentials" L. Larini, L. Lu, and G. A. Voth, Journal of Chemical Physics, 132, 164107 (2010). doi:10.1063/1.3394863

* "Fitting coarse-grained distribution functions through an iterative force-matching method" L. Lu, J. F. Dama, and G. A. Voth, Journal of Chemical Physics, 139, 21906 (2013). doi:10.1063/1.4811667

* "The Theory of Ultra-Coarse-Graining. 2. Numerical Implementation" A. Davtyan, J. F. Dama, A. V. Sinitskiy, and G. A. Voth, Journal of Chemical Theory and Computation, 10 (12), 5265–5275 (2014). doi:10.1021/ct500834t

* "The multiscale coarse-graining method. XI. Accurate interactions based on the centers of charge of coarse-grained sites", Z. Cao and G. A. Voth, Journal of Chemical Physics, 143, 243116 (2015). doi:10.1063/1.4933249

* "Predicting the Sensitivity of Multiscale Coarse-Grained Models to their Underlying Fine-Grained Model Parameters" J. W. Wagner, J. F. Dama, and G. A. Voth, Journal of Chemical Theory and Computation, 11(8), 3547–3560 (2015). doi:10.1021/acs.jctc.5b00180

* "A Direct Method for Incorporating Experimental Data into Multiscale Coarse-Grained Models", T. Dannenhoffer-Lafage, A. D. White, and G. A. Voth, Journal of Chemical Theory and Computation, 12(5), 2144–2153 (2016). doi:10.1021/acs.jctc.6b00043

* "The Theory of Ultra-Coarse-Graining. 3. Coarse-Grained Sites with Rapid Local Equilibrium of Internal States", J. F. Dama, J. Jin, and G. A. Voth, Journal of Chemical Theory and Computation, 2017 ASAP Article. doi: 10.1021/acs.jctc.6b01081

* "Extending the Range and Physical Accuracy of Coarse-grained Models: Order Parameter Dependent Interactions", J. W. Wagner, T. Dannenhoffer-Lafage, J. Jin, and G. A. Voth, Journal of Chemical Physics, 147, 044113 (2017). doit:10.1063/1.4995946
