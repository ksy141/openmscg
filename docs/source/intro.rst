Introduction
============

The open multiscale coarse-graining (openMSCG) python library provides various
systematic, bottom-up methods to calculate an effective CG force-field based on
rigorous statistical mechanics. This includes approximating the many-body
potential of mean force by minimizing the mean squared error between CG forces
and the mapped fine-grained forces (a.k.a, "force matching"). 

**This package also provides support for the following methods:**

* Relative entropy

* Iterative boltzmann inversion

* Rapid local equilibrium Ultra Coarse-Graining

**The terminology and algorithms used in the code are explained the publications
for these methods, which include:**

* "Efficient, Regularized, and Scalable Algorithms for Multiscale Coarse-Graining", L. Lu, S. Izvekov, A. Das, H. C. Andersen, and G. A. Voth, Journal of Chemical Theory and Computation, 6(3), 954-965 (2010). doi:10.1021/ct900643r

* "Multiscale coarse graining of liquid-state systems", S. Izvekov, and G. A. Voth, The Journal of Chemical Physics, 123, 134105 (2005). doi:10.1063/1.2038787

* "The multiscale coarse-graining method. II. Numerical implementation for coarse-grained molecular models", W. G. Noid, P. Liu, Y. Wang, J.-W. Chu, G. S. Ayton, S. Izvekov, H. C. Andersen, and G. A. Voth, The Journal of Chemical Physics, 128, 244115 (2008). doi:10.1063/1.2938857

* "The theory of ultra-coarse-graining. 2. Numerical implementation", A. Davtyan, J. F. Dama, A. V. Sinitskiy, and G. A. Voth, The Journal of Chemical Theory and Computation, 10, 5265–5275 (2014). doi:10.1021/ct500834t

* "The theory of ultra-coarse-graining. 3. Coarse-grained sites with rapid local equilibrium of internal states", J. F. Dama, J. Jin, and G. A. Voth, The Journal of Chemical Theory and Computation, 13, 1010-1022 (2017). doi: 10.1021/acs.jctc.6b01081

* "Coarse-graining errors and numerical optimization using a relative entropy framework", A. Chaimovich and M. Shell, The Journal of Chemical Physics, 134, 094112 (2011). doi:10.1063/1.3557038


Funding
#######

This project is funded by NSF/SSE (No.5-27425), and focuses on the sustainable
development of a contemporary CG/UCG molecular simulation platform with broad
applicability to the fields of chemistry, biology, and materials science
in accordance with the National Strategic Supercomputing Initiative (e.g.
Principle 1, Objective 2), and the Vision and Strategy for Software for
Science, Engineering, and Education (NSF 12-113, Strategies 1-4).

