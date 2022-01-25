function [syssch,Q,Z,ev] = slgsrsf(sys)
%SLGSRSF Transform the pair (A,E) of a descriptor system to a
%        real generalized Schur form.
%
%        [SYSSCH,Q,Z,EV] = SLGSRSF(SYS)  applies an orthogonal
%        equivalence transformation to the pair (A,E) of a descriptor
%        system, SYS = (A,E,B,C), so that the transformed pair (in SYSSCH),
%        (A,E) <- (Q'*A*Z,Q'*E*Z), is in a real generalized Schur form.
%        The transformation is also applied to the matrices B and C.
%        The poles of the descriptor system are returned in EV.
%
%        SYSSCH = SLGSRSF(SYS)  performs the same computations,
%        but Q, Z, and EV are not returned, although they are computed.
%
%        Note: This function cannot work with a singular matrix E. 
%              Use GSYSTRA instead.
%
%        See also SLGSBAL, SLGSHES, SLGSQRQ, SLGSSVD, GSYSTRA
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 30-04-2003, 03-03-2009.
%

[A,B,C,D,E] = dssdata( sys );
if norm( E - eye( size( E ) ), 1 ) == 0,
   % Standard system: use just real Schur form.
   if nargout <= 1
      [A,B,C] = systra( 2,A,B,C );
   else  
      [A,B,C,ev,Z] = systra( 2,A,B,C );  Q = Z';
      if isempty( ev ),  ev = [];  end
   end
   syssch = ss( A,B,C,D,sys );
%
else
   % Descriptor system: use real generalized Schur form.
   if nargout <= 1
      [A,E,B,C] = gsystra( 6,A,E,B,C );
   else  
      [A,E,B,C,Q,Z,ev,db] = gsystra( 6,A,E,B,C );
      if nargout == 4,  
         ev = ev./db;  if isempty( ev ),  ev = [];  end
      end
   end
   syssch = dss( A,B,C,D,E,sys );
end
%
% end slgsrsf
