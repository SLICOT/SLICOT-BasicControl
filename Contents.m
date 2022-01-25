% SLBASIC - SLICOT Basic Systems and Control Toolbox.
% RELEASE 2.0 of SLICOT Toolboxes.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% This toolbox includes M- and MEX-functions for basic
% computations in systems and control.
%
% MATLAB M-files (except help files).
% -----------------------------------
%
% SLICOT-based MATLAB functions.
%
%  Basic tools for standard and generalized state space systems and
%  transfer matrix factorizations
% 
%   a. Canonical forms and system transformations.
%   slconf      - Controllability staircase form of a system SYS = (A,B,C). 
%   slobsf      - Observability staircase form of a system SYS = (A,B,C).
%   slminr      - Minimal realization of a system SYS = (A,B,C).
%   ckstair     - Check that a system is in a staircase form.
%   slsbal      - Balance the system matrix for a state-space system SYS = (A,B,C). 
%   slsdec      - Additive spectral decomposition of a system SYS = (A,B,C) with
%                 respect to a given stability domain. 
%   lcf         - Left coprime factorization with prescribed stability degree.
%   lcfid       - Left coprime factorization with inner denominator.
%   rcf         - Right coprime factorization with prescribed stability degree.
%   rcfid       - Right coprime factorization with inner denominator.
%   slsorsf     - Transform the state matrix of a state space system
%                 to a specified eigenvalue-ordered real Schur form. 
%   slsrsf      - Transform the state matrix A to a real Schur form. 
%   slc2d       - Bilinear transformation of a continuous-time system 
%                 to a discrete-time system.
%   sld2c       - Bilinear transformation of a discrete-time system 
%                 to a continuous-time system.
%   cf2ss       - Construct the state-space representation of a system from
%                 the factors of its left or right coprime factorization.
%   polysf      - Compute the spectral factorization of a real polynomial.
%
%   b. System inter-connections. 
%   slapp       - Append two systems in state-space form.
%   slfeed      - Feedback inter-connection of two systems in
%                 state-space form.
%   slpar       - Parallel inter-connection of two systems in
%                 state-space form.
%   slser       - Series inter-connection of two systems in
%                 state-space form.
%   slspar      - Rowwise inter-connection of two systems in
%                 state-space form.
%   slosfeed    - Closed-loop system for a mixed output and state
%                 feedback control law.
%   slofeed     - Closed-loop system for an output feedback
%                 control law.
%
%   c. Dual and inverse.
%   sldual      - Dual of a standard system.
%   slinv       - Inverse of a standard system.
%
%   d. System norms. 
%   slH2norm    - H2/L2 norm of a system.
%   slHknorm    - Hankel-norm of a stable projection of a system.
%   slinorm     - L-infinity system norm.
%   slstabr     - complex stability radius.
%
%   e. Canonical forms and system transformations for descriptor systems.
%   slgconf     - Reduce a descriptor system to controllable staircase form.
%   slgobsf     - Reduce a descriptor system to observable staircase form.
%   slgminr     - Reduce a descriptor system to an irreducible form.
%   slgsbal     - Balance the system matrix for a descriptor system.
%   slgsHes     - Transform the pair (A,E) of a descriptor system to a
%                 generalized Hessenberg form.
%   slgsQRQ     - Transform the pair (A,E) of a descriptor system to a
%                 QR- or RQ-coordinate form.
%   slgsrsf     - Transform the pair (A,E) of a descriptor system to a
%                 real generalized Schur form.
%   slgsSVD     - Transform the pair (A,E) of a descriptor system to a
%                 singular value decomposition (SVD) or SVD-like coordinate
%                 form.
%   nrank       - Compute the normal rank of the transfer-function matrix
%                 of a standard system.
%   polzer      - Compute the normal rank, poles, zeros, and the Kronecker
%                 structure of the system pencil for a standard or descriptor
%                 system.
%   slpole      - Compute the poles of a standard or descriptor system.
%   slzero      - Compute the normal rank, zeros, and the Kronecker
%                 structure of the system pencil for a standard or descriptor
%                 system.
%
%   f. Matrix exponential and integral.
%   slexpe      - Compute the matrix exponential basically using an
%                 eigenvalue/eigenvector decomposition technique, but
%                 a diagonal Pade approximant with scaling and squaring,
%                 if the matrix appears to be defective.
%   slexpm      - Compute the matrix exponential using a diagonal Pade
%                 approximant with scaling and squaring.
%   slexpi      - Compute the matrix exponential and optionally its 
%                 integral, using a Pade approximation of the integral.
%
%   g. Pole assignment.
%   deadbt      - Construct the minimum norm feedback matrix performing 
%                 "deadbeat control" on an (A,B)-pair.
%   pass        - Perform (partial) pole assignment.
%
%   h. Kalman filter.
%   convKf      - Compute one recursion of the conventional Kalman filter
%                 equations.
%   srcf        - Compute a combined measurement and time update of one
%                 iteration of the time-varying or time-invariant Kalman
%                 filter in the Square Root Covariance Form.
%   srif        - Compute a combined measurement and time update of one
%                 iteration of the time-varying or time-invariant Kalman
%                 filter in the Square Root Information Form.
%   tLSfilt     - Compute a combined measurement and time update of
%                 the recursive least-squares filter. It is intended
%                 for testing the gateway Kfiltupd for task = 3.
%
%   i. Riccati equations.
%   slcaregs    - Solve CARE with generalized Schur method on an extended
%                 matrix pencil.
%   slcares     - Solve CARE with Schur method.
%   slcaresc    - Solve CARE with refined Schur method and estimate condition. 
%   sldaregs    - Solve DARE with generalized Schur method on an extended
%                 matrix pencil. 
%   sldares     - Solve DARE with Schur method. 
%   sldaresc    - Solve DARE with refined Schur method and estimate condition. 
%   sldaregsv   - Solve DARE with generalized Schur method on a symplectic 
%                 pencil. 
%   slgcare     - Solve descriptor CARE with generalized Schur method. 
%   slgdare     - Solve descriptor DARE with generalized Schur method. 
%   carecond    - Estimate reciprocal condition number of a CARE.
%   darecond    - Estimate reciprocal condition number of a DARE.
%
%   j. Sylvester and Lyapunov-like equations.
%   slsylv      - Solve continuous-time Sylvester equations. 
%   sldsyl      - Solve discrete-time Sylvester equations. 
%   sllyap      - Solve continuous-time Lyapunov equations. 
%   slstei      - Solve Stein equations. 
%   slstly      - Solve stable continuous-time Lyapunov equations. 
%   slstst      - Solve stable Stein equations. 
%   slgesg      - Solve generalized linear matrix equation pairs. 
%   slgely      - Solve generalized continuous-time Lyapunov equations. 
%   slgest      - Solve generalized Stein equations. 
%   slgsly      - Solve stable generalized continuous-time Lyapunov equations. 
%   slgsst      - Solve stable generalized Stein equations. 
%   lyapcond    - Estimate reciprocal condition number of a Lyapunov equation.
%   steicond    - Estimate reciprocal condition number of a Stein equation.
%
%   k. Output response.
%   dsimt       - Output response of a linear discrete-time system.
%                 The trajectories are stored column-wise.
%
%   l. Special numerical linear algebra computations.
%   bdiag       - Block diagonalization of a general matrix or a matrix
%                 in real Schur form.
%   persch      - Compute the periodic Hessenberg or periodic Schur decomposition
%                 of a matrix product.
%   Hameig      - Compute the eigenvalues of a Hamiltonian matrix.
%   habalance   - Symplectic scaling of a Hamiltonian matrix to improve
%                 eigenvalue accuracy.
%   haconv      - Convertions between storage representations for a
%                 Hamiltonian matrix.
%   haeig       - Eigenvalues of a Hamiltonian matrix using HAPACK approach.
%   hapvl       - Paige-Van Loan's form of a Hamiltonian matrix.
%   haschord    - Reordering the Schur form of a Hamiltonian matrix.
%   hastab      - Complete stable/unstable invariant subspace of a
%                 Hamiltonian matrix.
%   hasub       - Selected stable/unstable invariant subspace of a
%                 Hamiltonian matrix.
%   haurv       - Symplectic URV form of a general 2n-by-2n matrix.
%   haurvps     - Symplectic URV/periodic Schur form of a Hamiltonian matrix.
%   shbalance   - Symplectic scaling of a skew-Hamiltonian matrix to improve
%                 eigenvalue accuracy.
%   shconv      - Conversions between storage representations for a
%                 skew-Hamiltonian matrix.
%   Hesscond    - Compute an estimate of the reciprocal of the condition
%                 number of an upper Hessenberg matrix.
%   shheig      - Eigenvalues of a skew-Hamiltonian/(skew-)Hamiltonian matrix pencil.
%   shhstab     - Complete stable right deflating subspace of a skew-Hamiltonian/
%                 Hamiltonian matrix pencil. The stable companion subspace can
%                 also be returned for a pencil with the skew-Hamiltonian matrix 
%                 in factored form.
%   Hessl       - Solve a set of systems of linear equations with an
%                 upper Hessenberg coefficient matrix.
%   TLS         - Solve the Total Least Squares (TLS) problem using a 
%                 singular value decomposition (SVD) approach or a Partial
%                 SVD (PSVD) approach.
%
%  Benchmarks
% 
%   aredata     - Generate benchmark examples for algebraic Riccati equations.
%   ctdsx       - Generate benchmark examples for time-invariant, continuous-time,
%                 dynamical systems. 
%   ctlex       - Generate benchmark examples of (generalized) continuous-time 
%                 Lyapunov equations.
%   dtdsx       - Generate benchmark examples for time-invariant, discrete-time, 
%                 dynamical systems.
%   dtlex       - Generate benchmark examples of (generalized) discrete-time 
%                 Lyapunov equations.
%
%  Basic tools for structured matrix factorizations
% 
%   fstchol     - Factor a symmetric positive definite (block) Toeplitz matrix BT
%                 and solve associated linear systems, given the first (block) 
%                 row / column T of BT.
%   fstgen      - Factor a symmetric positive definite (block) Toeplitz matrix BT,
%                 computes the generator of inv(BT), and/or solve associated linear  
%                 systems using the Cholesky factor of inv(BT), given the first 
%                 (block) row / column T of BT.
%   fstsol      - Solve linear systems X*BT = B / BT*X = B, where BT is a symmetric 
%                 positive definite (block) Toeplitz matrix, given the first 
%                 (block) row / column T of BT. 
%   fstupd      - Factor and/or update a factorization of a symmetric positive 
%                 definite (block) Toeplitz matrix BT or BTA and solve associated
%                 linear systems, given the first (block) row / column T of BT or
%                 [ T Ta ] / [ T  ] of BTA.
%                            [ Ta ]
%   fstqr       - Compute the orthogonal-triangular decomposition of a (block)
%                 Toeplitz matrix BT and solve associated linear least-squares
%                 problems, given the first (block) column and the first (block)
%                 row of BT. It is assumed that the first MIN(SIZE(BT)) columns
%                 of BT have full rank.
%   fstlsq      - Solve linear least-squares problems min(B-BT*X) or find the
%                 minimum norm solution of BT'*Y = C where BT is a (block)
%                 Toeplitz matrix with full column rank, given the first (block)
%                 column and the first (block) row of BT.
%   fstmul      - Compute the matrix-vector products X = BT*B for a block Toeplitz
%                 matrix BT, given the first block column and the first block
%                 row of BT.
%
%  Data analysis
% 
%   sincos      - Sine or cosine transform of a real signal.
%   slDFT       - Discrete Fourier transform of a signal.
%   sliDFT      - Inverse discrete Fourier transform of a signal.
%   slHart      - Discrete Hartley transform of a real signal.
%   slconv      - Convolution of two real signals using either FFT or 
%                 Hartley transform.
%   sldeconv    - Deconvolution of two real signals using either FFT or 
%                 Hartley transform.
%   sncs        - Compute the sine or cosine transform of a real signal.
%   slwindow    - Anti-aliasing window applied to a real signal.
%
%  Mu optimal and H_infinity synthesis 
% 
%   slihinf     - Compute the H_infinity optimal controller for a state
%                 space model.
%   slimju      - Compute the mu optimal controller for a state space
%                 model, and the mu norm of the closed loop system.
%
% SLICOT-based MATLAB MEX-files.
% ------------------------------
%
%  Basic tools for standard and generalized state space systems and
%  transfer matrix factorizations.
% 
%   a. Canonical forms and system transformations.
%   syscom      - Compute controllability/observability forms, and minimal
%                 realization of a given system, based on staircase form
%                 reduction.
%   systra      - Compute system similarity transformations with scaling, block
%                 diagonal decomposition, or Schur form reduction of the state
%                 matrix.
%   condis      - Perform a transformation on the parameters (A,B,C,D) of a
%                 system, which is equivalent to a bilinear transformation
%                 of the corresponding transfer function matrix.
%   cfsys       - Construct the state-space representation of a system from
%                 the factors of its left or right coprime factorization.
%   specfact    - Compute the spectral factorization of a real polynomial.
%
%   b. System inter-connections.
%   sysconn     - Compute a state-space model (A,B,C,D) for various
%                 inter-connections of two systems given in state-space
%                 form.
%   sysfconn    - Compute, for a given state-space system (A,B,C,D), the
%                 closed-loop system (Ac,Bc,Cc,Dc) corresponding to the output,
%                 or mixed output and state, feedback control law.
%
%   c. Dual and inverse.
%   invert      - Compute the dual or inverse of a linear (descriptor)
%                 system.
%
%   d. System norms. 
%   Hnorm       - Compute various system norms (Hankel norm, H2 norm) and
%                 complex stability radius.
%   linorm      - Compute the L-infinity norm for a continuous-time or
%                 discrete-time system, possibly unstable and/or of the
%                 descriptor type.
%
%   e. Canonical forms and system transformations for descriptor systems.
%   gsyscom     - Transform a descriptor system, by equivalence
%                 transformations, to a controllable or observable staircase
%                 form, or to a reduced (controllable, observable, or
%                 irreducible) form.
%   gsystra     - Perform various equivalence transformations for descriptor
%                 systems with scaling, generalized Schur form, etc.
%   polezero    - Compute the normal rank, poles, zeros, and the Kronecker
%                 structure of the system pencil for a standard or descriptor
%                 system.
%   polezeroz   - Compute the normal rank, poles, zeros, and the Kronecker
%                 structure of the system pencil for a standard or descriptor
%                 complex system.
%   syscf       - Compute a left coprime factorization (LCF) or a right coprime 
%                 factorization (RCF) of a transfer-function matrix corresponding
%                 to a state-space system (A,B,C,D).
%   isprpr      - Check out the properness of the transfer function of
%                 a descriptor system.
%
%   f. Matrix exponential and integral.
%   slmexp      - Compute the matrix exponential and optionally its
%                 integral.
%
%   g. Pole assignment.
%   deadbeat    - Construct the minimum norm feedback matrix performing 
%                 "deadbeat control" on an (A,B)-pair.
%   polass      - Perform (partial) pole assignment.
%
%   h. Kalman filter.
%   Kfiltupd    - Compute a combined measurement and time update of one
%                 iteration of the Kalman filter.
%
%   i. Riccati equations.
%   aresol      - Solve algebraic Riccati equations.
%   aresolc     - Solve algebraic Riccati equations with condition and forward
%                 error estimates (refined Schur vector method).
%   garesol     - Solve descriptor algebraic Riccati equations.
%   arecond     - Estimate condition and forward errors for Lyapunov and 
%                 algebraic Riccati equations.
%
%   j. Sylvester and Lyapunov-like equations.
%   linmeq      - Solve linear matrix equations (Sylvester and Lyapunov).
%   genleq      - Solve generalized linear matrix equations.
%
%   k. Output response.
%   ldsimt      - Compute the output response of a linear discrete-time
%                 system.  The input and output trajectories are stored 
%                 column-wise (each column contains all inputs or outputs 
%                 measured at a certain time instant).
%
%   l. Special numerical linear algebra computations.
%   Hessol      - Analysis and solution of a system of linear equations with
%                 an upper Hessenberg coefficient matrix.
%   bldiag      - Perform block-diagonalization of a general matrix or a
%                 matrix in real Schur form.
%   perschur    - Compute the periodic Hessenberg or periodic Schur decomposition
%                 of a matrix product.
%   Hamileig    - Compute the eigenvalues of a Hamiltonian matrix using the
%                 square-reduced approach.
%   hapack_haeig- Compute the eigenvalues of a Hamiltonian matrix using the
%                 HAPACK approach.
%   HaeigZ      - Compute the eigenvalues of a complex Hamiltonian matrix.
%   skewHamil2eig   - Compute the eigenvalues of a real skew-Hamiltonian/
%                     skew-Hamiltonian pencil.
%   skewHamil2feig  - Compute the eigenvalues of a real skew-Hamiltonian/
%                     skew-Hamiltonian pencil (factored version).
%   skewHamildefl   - Compute the eigenvalues of a skew-Hamiltonian/Hamiltonian
%                     pencil and the right deflating subspace corresponding to the
%                     eigenvalues with strictly negative real part.
%   skewHamildeflf  - Compute the eigenvalues of a real skew-Hamiltonian/
%                     Hamiltonian pencil and the right deflating subspace
%                     corresponding to the eigenvalues with strictly negative
%                     real part (factored version).
%   skewHamildeflfZ - Compute the eigenvalues of a complex skew-Hamiltonian/
%                     Hamiltonian pencil and the right deflating subspace
%                     corresponding to the eigenvalues with strictly negative
%                     real part (factored version).
%   skewHamildeflZ  - Compute the eigenvalues of a complex skew-Hamiltonian/
%                     Hamiltonian pencil and the right deflating subspace
%                     corresponding to the eigenvalues with strictly negative
%                     real part.
%   skewHamileig    - Compute the eigenvalues and orthogonal decomposition of a
%                     skew-Hamiltonian/Hamiltonian pencil.
%   skewHamileigZ   - Compute the eigenvalues of a complex skew-Hamiltonian/
%                     Hamiltonian pencil.
%   symplURV    - Compute the eigenvalues and generalized symplectic URV
%                 decomposition of a skew-Hamiltonian/Hamiltonian pencil in
%                 factored form.
%   symplURVZ   - Compute the eigenvalues and generalized symplectic URV
%                 decomposition of a complex skew-Hamiltonian/Hamiltonian
%                 pencil in factored form.
%   TotalLS     - Solve the Total Least Squares (TLS) problem using a 
%                 singular value decomposition (SVD) approach or a Partial
%                 SVD (PSVD) approach.
%
%  Benchmarks
%   arebench    - Generate benchmark examples for algebraic Riccati equations.
% 
%  Basic tools for structured matrix factorizations.
%
%   fstoep      - Factor symmetric positive definite (block) Toeplitz matrices
%                 and/or solve associated linear systems.
%   fstoeq      - Compute orthogonal-triangular decompositions of (block)
%                 Toeplitz matrices and/or solve associated linear
%                 least-squares systems.
%
%  Data analysis.
% 
%   datana      - Perform various transforms of real or complex vectors.
%   
%  Mu optimal and H_infinity synthesis.
% 
%   conhin      - Compute an H-infinity or H2 controller for a continuous-time
%                 system.
%   dishin      - Compute an H-infinity or H2 controller for a discrete-time
%                 system.
%   clsdp       - Compute the controller for the Loop Shaping Design of a
%                 continuous-time system.
%   dlsdp       - Compute the controller for the Loop Shaping Design of a
%                 discrete-time system.
%   mucomp      - Compute the structured singular value.
%   dlsdz       - Loop Shaping Design of discrete-time systems.
%   muHopt      - Mu optimal or H_infinity controller design.
%
% Demonstrations.
% ---------------
%
%   basicdemo   - Basic SLICOT demonstration.
%   fstdemo     - Structured matrix decompositions demonstration.
%
% Note:
% -----
%
%   Command-line help functions are available for all MATLAB M- and 
%   MEX-functions included in this toolbox.
%

%  ..CONTRIBUTOR..
%
%   V. Sima, Katholieke Univ. Leuven, Belgium, Sept. 2001.
%
%   Revisions:
%   V. Sima,  April 25, 2005, March 9, 2009, Jan. 2022.

