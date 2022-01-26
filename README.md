# **SLICOT Basic Systems and Control Toolbox**  

## About 

The `SLICOT Basic Systems and Control Toolbox` (`SLICOT-BasicControl`) includes [SLICOT](http://slicot.org/)-based MATLAB and Fortran tools for solving efficiently and reliably various basic computational problems for linear time-invariant multivariable systems analysis and synthesis. Standard and generalised (descriptor) state space systems are covered.

The main functionalities of the toolbox include:

  * similarity and equivalence transformations for standard and descriptor systems
    - essential computations with structured matrices, including
    - eigenvalues of a Hamiltonian matrix or of a skew-Hamiltonian/Hamiltonian pencil
    - Periodic Hessenberg and periodic Schur decompositions
  * computations with (block) Toeplitz matrices and systems
  * analysis of standard and descriptor systems
  * solution of Lyapunov and Riccati equations with condition estimation
  * coprime factorization and spectral decomposition of transfer-function matrices

The toolbox main features are:

  *  computational reliability using square-root and balancing-free accuracy enhancing
  *   high numerical efficiency, using latest algorithmic developments, structure exploiting algorithms, and dedicated linear algebra tools
  *   flexibility and easy-of-use
  *   enhanced functionality, e.g, for controller reduction
  *   standardized interfaces

The programs have been extensively tested on various test examples and are fully documented.

## Requirements

The codes have been tested with MATLAB 2015b through 2021b. To use the functions, the Control System Toolbox must be installed in MATLAB running under 64-bit Windows 7, 8, 8.1 or 10. 

## License

* See [`LICENSE`](https://github.com/SLICOT/SLICOT-BasicControl/blob/master/LICENSE) for licensing information.

## References

Please cite `SLICOT-BasicControl` using at least one of the following references: 

* P. Benner, D. Kressner, V, Sima, and A. Varga, [The SLICOT Toolboxes - a Survey](http://slicot.org/objects/software/reports/SLWN2009-1.pdf), _SLICOT Working Note 2009-1, August 2009._
* P. Benner, D. Kressner, V. Sima, A. Varga, Die SLICOT-Toolboxen für Matlab - The SLICOT Toolboxes for Matlab (in German), _at – Automatisierungstechnik, 58 (2010)._

