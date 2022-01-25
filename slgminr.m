function [sysr,infred,s] = slgminr(sys,job,tol,systyp,bal)
%SLGMINR Transform a descriptor system (A-lambda*E,B,C) to a
%        controllable, observable, or irreducible form.
%
%        SYSR = SLGMINR(SYS)  computes a reduced (controllable,
%        observable, or irreducible) descriptor representation
%        (Ar-lambda*Er,Br,Cr) for an original descriptor representation
%        (A-lambda*E,B,C). The pencil Ar-lambda*Er is in an upper block
%        Hessenberg form, with either Ar or Er upper triangular. The
%        transformed system is returned in the LTI object SYSR. 
%
%        [SYSR,INFRED,S] = SLGMINR(SYS,JOB,TOL,SYSTYP,BAL)  has
%        additional input and output arguments.
%
%        JOB is an integer scalar specifying what to do, as follows: 
%        JOB = 0 :  remove both the uncontrollable and unobservable parts
%                   to get an irreducible descriptor representation;
%        JOB = 1 :  remove the uncontrollable part only to get a
%                   controllable descriptor representation;
%        JOB = 2 :  remove the unobservable part only to get an
%                   observable descriptor representation.
%        Default:  JOB = 0.
%
%        TOL is the tolerance to be used in rank determinations when
%        transforming (A-lambda*E,B,C). If TOL > 0, then the given value
%        of TOL is used as a lower bound for reciprocal condition numbers;
%        a (sub)matrix whose estimated condition number is less than
%        1/tol is considered to be of full rank.
%        If TOL <= 0, the default tolerance toldef = eps*prod(size(A))
%        is used instead, where eps is the machine precision. TOL < 1.
%        Default :  TOL = 0.
%
%        SYSTYP is an integer scalar indicating the type of descriptor
%        system algorithm to be applied according to the assumed
%        transfer-function matrix as follows:
%        SYSTYP = 0 :  rational transfer-function matrix;
%        SYSTYP = 1 :  proper (standard) transfer-function
%                      matrix;
%        SYSTYP = 2 :  polynomial transfer-function matrix.
%        Default:  SYSTYP = 0.
%
%        BAL is an integer scalar specifying whether the user wishes to
%        preliminarily scale the system (A-lambda*E,B,C) as follows:
%        BAL = 0 :  perform scaling;
%        BAL = 1 :  do not perform scaling.
%        Default:  BAL = 0.
%
%        The matrix Ar of order nr is upper triangular if SYSTYP = 0 or 2.
%        If SYSTYP = 1 and JOB = 1, the matrix [Br Ar] is in a
%        controllable staircase form.
%        If SYSTYP = 1 and JOB = 0 or 2, the matrix ( Ar ) is in an
%                                                   ( Cr )
%        observable staircase form.
%        The block structure of staircase forms is contained
%        in the leading INFRED(7) elements of the vector S.
%
%        The matrix Er has INFRED(6) nonzero sub-diagonals.
%        If at least for one k = 1,...,4, INFRED(k) >= 0, then the
%        resulting Er is structured being either upper triangular
%        or block Hesseneberg, in accordance to the last
%        performed order reduction phase (see Method).
%        The block structure of staircase forms is contained
%        in the leading INFRED(7) elements of the vector S.
%
%        If JOB = 1, only the first S(1) rows of B are nonzero.
%
%        If JOB = 0, or JOB = 2, only the last S(1) columns
%        (in the first nr columns) of C are nonzero.
%
%        INFRED is an integer array of dimension 7, containing
%        information on performed reduction and on structure of
%        resulting system matrices as follows:
%        INFRED(k) >= 0 (k = 1, 2, 3, or 4) if Phase k of reduction
%                       (see Method) has been performed. In this
%                       case, INFRED(k) is the achieved order
%                       reduction in Phase k.
%        INFRED(k) < 0  (k = 1, 2, 3, or 4) if Phase k was not
%                       performed.
%        INFRED(5)  -   the number of nonzero sub-diagonals of A.
%        INFRED(6)  -   the number of nonzero sub-diagonals of E.
%        INFRED(7)  -   the number of blocks in the resulting
%                       staircase form at last performed reduction
%                       phase. The block dimensions are contained
%                       in the first INFRED(7) elements of S.
%
%        S is an integer vector of dimension INFRED(7) containing 
%        the orders of the diagonal blocks of Ar-lambda*Er.
%
%        For standard systems (i.e., with identity E), the output 
%        arguments INFRED and S are not set.
%
%        Method
%        The order reduction is performed in 4 phases:
%        Phase 1: Eliminate all finite uncontrolable eigenvalues.
%        The resulting matrix ( Br Ar ) is in a controllable
%        staircase form, and Er is upper triangular.
%        This phase is performed if JOB = 0 or 1 and SYSTYP = 0
%        or 1.
%        Phase 2: Eliminate all infinite and finite nonzero uncontrollable
%        eigenvalues. The resulting matrix ( Br Er ) is in a
%        controllable staircase form, and Ar is upper triangular.
%        This phase is performed if JOB = 0 or 1 and SYSTYP = 0
%        or 2.
%        Phase 3: Eliminate all finite unobservable eigenvalues.
%        The resulting matrix ( Ar ) is in an observable
%                             ( Cr )
%        staircase form, and Er is upper triangular.
%        This phase is performed if JOB = 0 or 2 and SYSTYP = 0
%        or 1.
%        Phase 4: Eliminate all infinite and finite nonzero unobservable
%        eigenvalues. The resulting matrix ( Er ) is in an
%                                          ( Cr )
%        observable staircase form, and Ar is upper triangular.
%        This phase is performed if JOB = 0 or 2 and SYSTYP = 0
%        or 2.
%
%        Note: This function cannot work with a singular matrix E.
%              Use GSYSCOM instead.
%
%        See also SLGCONF, SLGOBSF, GSYSCOM
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 15-05-2003.
%
%        Revisions: 03-03-2009.

ni = nargin;
%
if ni == 1,
    job    = 0;
    tol    = 0;
    systyp = 0;
    bal    = 0;
elseif ni == 2,
    tol    = 0;
    systyp = 0;
    bal    = 0;
elseif ni == 3,
    systyp = 0;
    bal    = 0;
elseif ni == 4,
    bal    = 0;
end
%
if job    < 0 || job    > 2,  error('Wrong value for JOB');     end
if systyp < 0 || systyp > 2,  error('Wrong value for SYSTYP');  end
%
[A,B,C,D,E] = dssdata( sys );
if norm( E - eye( size( E ) ), 1 ) == 0,
   % Standard system. The arguments INFRED and S are not set.
   sysr = slminr(sys,tol,bal);
   infred = [ -1; -1; -1; -1; 0; 0; 0 ];  s = [];
   %
else
   % Descriptor system.
   flag = [ job; systyp; bal; tol ];
   %
   [A,E,B,C,infred,s] = gsyscom( 3,A,E,B,C,flag );
   sysr = dss( A,B,C,D,E,sys );
end
%
% end slgminr
