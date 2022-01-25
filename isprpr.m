% ISPRPR.F - Gateway function to check whether the transfer function of
%            a descriptor system with given generalized state space
%            realization is proper or not, using the SLICOT Library
%            routine AB13ID.
%
%   [prop(,Aout,Eout,Bout,Cout,warn)] =
%         isprpr(A,E,B,C(,jobsys,jobeig,equil,cksing,restor,update,tol))
%   [prop(,warn)] =
%         isprpr(A,E,B,C(,jobsys,0,equil,cksing,restor,0,tol))
%
%   ISPRPR checks whether the transfer function
%
%                                   -1
%     G(lambda) := C*( lambda*E - A ) *B
%
%   of a given linear time-invariant descriptor system with generalized
%   state space realization (lambda*E-A,B,C) is proper.
%   Optionally, the system (lambda*E-A,B,C) is reduced to an equivalent
%   one (lambda*Er-Ar,Br,Cr) with only controllable and observable
%   eigenvalues in order to use it for a subsequent L_inf-norm
%   computation; or the system is reduced to an equivalent one
%   (lambda*Er-Ar,Br,Cr) without uncontrollable and unobservable
%   infinite eigenvalues. In this case, if update = 0, the system is
%   only checked for properness, the reduced system is not returned.
%
% Description of input parameters:
%   A       - the n-by-n system state matrix A.
%   E       - the n-by-n system descriptor matrix E.
%   B       - the n-by-m system input matrix B.
%   C       - the p-by-n system output matrix C.
%   jobsys  - (optional) specifies if the system is already in the
%             reduced form as stated in jobeig.
%             = 0 : system is not in reduced form;
%             = 1 : system is in reduced form.
%             Default: jobsys = 0.
%   jobeig  - (optional) specifies which kind of eigenvalues of the
%             system pencil lambda*E-A should be removed if jobsys = 0.
%             = 0 : remove only infinite uncontrollable and unobservable
%                   eigenvalues;
%             = 1 : remove all uncontrollable and unobservable
%                   eigenvalues.
%             Default: jobeig = 0.
%   equil   - (optional) specifies whether the user wishes to
%             preliminarily equilibrate the system (lambda*E-A,B,C):
%             = 0 : do not perform equilibration;
%             = 1 : perform equilibration (scaling).
%             Default: equil = 0.
%   cksing  - (optional) specifies whether the user wishes to check if
%             the pencil (A-lambda*E) is singular:
%             = 0 : do not check singularity;
%             = 1 : check singularity.
%             If the pencil is singular, the reduced system computed for
%             cksing = 0 can have completely different eigenvalues than
%             the given system.
%             The test is performed only if jobsys = 0.
%             Default: cksing = 1.
%   restor  - (optional) specifies whether the user wishes to save the
%             system matrices before each reduction phase if jobsys = 0,
%             and restore them if no order reduction took place:
%             = 0 : do not save the matrices;
%             = 1 : save and restore.
%             While saving and restoring can sometimes give more
%             accurate reduced systems (with poles/zeros closer to the
%             original ones), this is not always true. Use restor = 0
%             for testing the properness.
%             Default: restor = 0.
%   update  - (optional) specifies whether the user wishes to update the
%             matrices A, B, and C if jobeig = 0:
%             = 0 : do not update the matrices A, B and C;
%             = 1 : update the matrices A, B and C.
%             Default: update = 0.
%   tol     - (optional) vector of length 3 containing the tolerance
%             used to set the accuracy in determining ranks (in tol(1)),
%             the tolerance used for checking pencil singularity,
%             when cksing = 1, and/or the singularity of the matrices A
%             and E, when cksing = 0 (in tol(2)), and the threshold
%             value for magnitude of the matrix elements, if equil = 1:
%             elements with magnitude less than or equal to tol(3) are
%             ignored for scaling.
%             Default: tol(1) = n*n*(epsilon_machine),
%                      tol(2) =  10*(epsilon_machine),
%                      tol(3) = max(norm_1(A,E,B,C))*(epsilon_machine),
%             where epsilon_machine is the relative machine precision.
%
% Description of output parameters:
%   prop    - indicates whether the transfer function of the system is
%             proper or not.
%             = 0 : the transfer function is improper;
%             = 1 : the transfer function is proper.
%   Aout    - (optional) if jobeig = 1 or update = 1, the nr-by-nr
%             reduced or transformed system state matrix Aout.
%   Eout    - (optional) if jobeig = 1 or update = 1, the nr-by-nr
%             reduced or transformed system descriptor matrix Eout.
%   Bout    - (optional) if jobeig = 1 or update = 1, the nr-by-m
%             reduced or transformed system input matrix Bout.
%   Cout    - (optional) if jobeig = 1 or update = 1, the p-by-nr
%             reduced or transformed system output matrix Cout.
%   warn    - (optional) warning indicator which indicates the quality
%             of rank decisions
%             = 0 : rank decisions can be assumed to be correct;
%             = 1 : the tolerance tol(1) is very close to a singular
%                   value of a matrix whose rank should be computed;
%                   rank decisions might be incorrect.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2012.
%
% Revisions:
%   -
%
