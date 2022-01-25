% SYSFCONN.F - MEX-function to compute, for a given state-space
%              system (A,B,C,D), the closed-loop system (Ac,Bc,Cc,Dc)
%              corresponding to the output, or mixed output and state,
%              feedback control law, using SLICOT routines AB05RD and
%              AB05SD.
%
%   [Ao,Bo,Co(,Do,rcondi)] = SYSFCONN(task,jobd,fbtype,alpha(,beta)
%                                     A,B,C(,D)(,G,H)(,F)(,K))
%
%   [Ao,Bo,Co(,Do,rcondi)] = SYSFCONN(1,jobd,fbtype,alpha,beta,A,B,C(,D)
%                                     G,H(,F)(,K))
%   [Ao,Bo,Co(,Do,rcondi)] = SYSFCONN(2,jobd,fbtype,alpha,A,B,C(,D)(,F))
%
%   SYSFCONN computes a state-space model (A,B,C,D) for the various
%   feedback inter-connections of two systems, according to the value 
%   of task:
%
%   task = 1: To construct, for a given state space system (A,B,C,D),
%   the closed-loop system (Ac,Bc,Cc,Dc) corresponding to the mixed
%   output and state feedback control law
%
%          u = alpha*F*y + beta*K*x + G*v,
%          z = H*y.                                                  (1)
%
%   task = 2: To construct, for a given state space system (A,B,C,D),
%   the closed-loop system (Ac,Bc,Cc,Dc) corresponding to the output
%   feedback control law
%
%          u = alpha*F*y + v.                                        (2)
%
%   Description of input parameters: 
%   task   - integer specifying the computations to be performed.
%            task = 1 :  compute the closed-loop system corresponding to
%                        a mixed output and state feedback control law;
%            task = 2 :  compute the closed-loop system corresponding to
%                        an output feedback control law.
%   jobd   - integer specifying whether or not a non-zero matrix D
%            appears in the given state space model:
%            jobd = 0 :  D is assumed a zero matrix;
%            jobd = 1 :  D is present.
%   fbtype - integer specifying the type of the feedback law as follows:
%            fbtype = 1 :  Unitary output feedback (F = I);
%            fbtype = 2 :  General output feedback.
%   alpha  - the real coefficient alpha in the output feedback law.
%   beta   - the real coefficient beta in the state feedback law.
%   A      - the n-by-n state dynamics matrix A.
%   B      - the n-by-m input/state matrix B.
%   C      - the p-by-n state/output matrix C.
%   D      - (optional) if jobd = 1, the p-by-m input/output matrix D.
%            If jobd = 0, this parameter must not be given.
%   G      - the m-by-mv system input scaling matrix G.
%   H      - the pz-by-p system output scaling matrix H.
%   F      - (optional) if fbtype = 2, the m-by-p output feedback
%            matrix F. If fbtype = 1, then the feedback matrix is
%            assumed to be an m-by-m identity matrix. In this case,
%            or if alpha = 0, the parameter F must not be given.
%   K      - (optional) if beta <> 0, the m-by-n state feedback
%            matrix K. If beta = 0, the parameter K must not be given.
%
%   Description of output parameters: 
%   Ao     - the n-by-n state dynamics matrix Ac of the closed-loop
%            system.
%   Bo     - the n-by-mb input/state matrix Bc of the closed-loop
%            system, where mb = mv if task = 1, and mb = m, is task = 2.
%   Co     - the pc-by-n output/state matrix Cc of the closed-loop
%            system, where pc = pz if task = 1, and pc = p, is task = 2.
%   Do     - (optional) if jobd = 1, the pc-by-mb input/output
%            matrix Dc of the closed-loop system.
%   rcondi - (optional) the reciprocal condition number of the matrix
%            I - alpha*D*F.
%
%   Method:
%   task = 1:  The matrices of the closed-loop system have the
%              expressions:
%
%     Ac = AC + beta*BC*K,      Bc = BC*G,
%     Cc = H*(CC + beta*DC*K),  Dc = H*DC*G,
%
%     where
%
%     AC = A + alpha*B*F*E*C,   BC = B + alpha*B*F*E*D,
%     CC = E*C,                 DC = E*D,
%
%     with E = (I - alpha*D*F)**-1.
%
%   task = 2:  The matrices of the closed-loop system have the
%              expressions:
%
%     Ac = A + alpha*B*F*E*C,  Bc = B + alpha*B*F*E*D,
%     Cc = E*C,                Dc = E*D,
%
%     where E = (I - alpha*D*F)**-1.
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
