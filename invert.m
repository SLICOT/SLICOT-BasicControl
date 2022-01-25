% INVERT.F   - MEX-function to compute the dual or inverse of a
%              linear (descriptor) system, using SLICOT routines
%              AB07MD, AB07ND, and AG07BD.
%
%   [Ao(,Eo),Bo,Co(,Do)(,rcondD)] = INVERT(task,A,B,C(,D)(,E))
%
%   [Ao,Bo,Co(,Do)]        = INVERT(1,A,B,C(,D))
%   [Ao,Bo,Co,Do(,rcondD)] = INVERT(2,A,B,C,D)
%   [Ao,Eo,Bo,Co,Do]       = INVERT(3,A,B,C,D(,E))
%
%   INVERT computes the dual of a standard system or the inverse of a
%   standard system or a descriptor system, according to the value of
%   task:
%
%   task = 1:  To compute the dual of a standard system in state-space
%   representation, i.e., if the m-input/p-output system is (A,B,C,D),
%   its dual is simply the p-input/m-output system (A',C',B',D').
%
%   task = 2:  To compute the inverse of a standard system.
%   The matrices of the inverse system are computed with the formulas
%                    -1              -1         -1           -1
%        Ai = A - B*D  *C,  Bi = -B*D  ,  Ci = D  *C,  Di = D  .     (1)
%
%   task = 3:  To compute the inverse (Ai-lambda*Ei,Bi,Ci,Di) of a given
%   descriptor system (A-lambda*E,B,C,D). The matrices of the inverse
%   system are computed with the formulas
%
%             ( E  0 )        ( A  B )         (  0 )
%        Ei = (      ) , Ai = (      ) ,  Bi = (    ),
%             ( 0  0 )        ( C  D )         ( -I )
%
%        Ci = ( 0  I ),  Di = 0.                                     (2)
%
%   Description of input parameters: 
%   task   - integer specifying the computations to be performed.
%            task = 1 :  compute the dual of a standard system;
%            task = 2 :  compute the inverse of a standard system;
%            task = 3 :  compute the inverse of a descriptor system.
%   A      - the n-by-n state dynamics matrix A.
%   B      - the n-by-m input/state matrix B.
%   C      - the p-by-n state/output matrix C. If task >= 2, p = m.
%   D      - the p-by-m input/output matrix D. If task >= 2, p = m.
%            If task = 1 and D is zero, this argument may be omitted.
%   E      - the n-by-n descriptor matrix E. When E is an identity
%            matrix, this argument may be omitted.
%
%   Description of output parameters: 
%   Ao     - If task = 1, the n-by-n dual state dynamics matrix A'.
%            If task = 2, the n-by-n state matrix Ai of the inverse
%            system (1).
%            If task = 3, the (n+m)-by-(n+m) state matrix Ai of the
%            inverse system (2).
%   Eo     - If task = 3 the (n+m)-by-(n+m) descriptor matrix Ei of the
%            inverse system (2).
%   Bo     - If task = 1, the n-by-p dual input/state matrix C'.
%            If task = 2, the n-by-m input matrix Bi of the inverse
%            system (1).
%            If task = 3, the (n+m)-by-m input matrix Bi of the
%            inverse system (2).
%   Co     - If task = 1, the m-by-n dual state/output matrix B'.
%            If task = 2, the m-by-n output matrix Ci of the inverse
%            system (1).
%            If task = 3, the m-by-(n+m) output matrix Ci of the
%            inverse system (2).
%   Do     - If task = 1, the m-by-p dual input/output matrix D'.
%            If task = 2, the m-by-m input/output matrix Di of the
%            inverse system (1).
%            If task = 3, the m-by-m input/output (zero) matrix Di of
%            the inverse system (2).
%   rcondD - (optional) the estimated reciprocal condition number of the
%            input/output matrix D of the original system.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, July 2003.
%
% Revisions:
%   -
%
