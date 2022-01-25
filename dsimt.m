function [Y,x] = dsimt(A,B,C,D,U,x0,HessA)
% DSIMT  Computes the output response of a linear discrete-time system.
%        The system state matrix may be given as an upper or lower
%        Hessenberg matrix.
%
%        [Y,x] = DSIMT(sys,U,x0)  computes the output vector sequence 
%        y(1), y(2),..., y(t) and the final state x of a discrete-time
%        state-space system, SYS = (A,B,C,D) (an ss object), given the
%        input vector sequence u(1), u(2),..., u(t) in U, and the initial
%        state vector x0. U and Y are matrices with t columns (t is the 
%        number of samples), and as many rows as inputs and outputs,
%        respectively.
%
%        [Y,x] = DSIMT(sys,U)  uses x0 = 0 as initial state.
%
%        [Y,x] = DSIMT(A,B,C,D,U,x0)  uses the system matrices instead sys.
%        [Y,x] = DSIMT(A,B,C,U,x0)    assumes that matrix D is zero.
%
%        [Y,x] = DSIMT(sys,U,x0,HessA)      enable to specify a Hessenberg
%        [Y,x] = DSIMT(A,B,C,D,U,x0,HessA)  form for the system state
%        [Y,x] = DSIMT(A,B,C,U,x0,HessA)    matrix.  
%
%        HessA is a scalar indicating whether the matrix A is general or
%        in an upper/lower Hessenberg form:
%        HessA = 0 :  general matrix;
%        HessA = 1 :  upper Hessenberg matrix;
%        HessA = 2 :  lower Hessenberg matrix.
%        Default: HessA = 0.
%
%        See also LDSIMT
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, April 2003.
%
%        Revisions:
%        V. Sima, March 2009.
%

ni = nargin;
if ni < 2,
   error( 'DSIMT needs at least 2 input parameters' )
end
if isa( A, 'lti' ),
   % Get the system matrices of the ss object, and the remaining parameters.
   % General call     [Y,x] = DSIMT(A,B,C,D,U,x0,HessA);
   % Special call     [Y,x] = DSIMT(sys,U,x0,HessA); 
   %
   if A.Ts == 0,
      error( 'The system SYS must be a discrete-time system' )
   end
   [ As, Bs, Cs, Ds ] = ssdata( A );
   [ l, m ] = size( A );  n = size( As, 1 );
   if ~( size( B, 1 ) == m ),
      error( 'The matrix U must have as many rows as inputs' )
   end
   if ni == 2,  
      % Special call     [Y,x] = DSIMT(sys,U);  Below, B is U ! 
      [ Y, x ] = ldsimt( As, Bs, Cs, Ds, B );
   elseif ni == 3,
      if numel( C ) < n,
         error( 'The initial state vector is too short' )
      end
      % Special call     [Y,x] = DSIMT(sys,U,x0);  Below, B is U and C is x0 ! 
      [ Y, x ] = ldsimt( As, Bs, Cs, Ds, B, C );
   elseif ni == 4,
      if numel( C ) < n,
         error( 'The initial state vector is too short' )
      end
      % Special call     [Y,x] = DSIMT(sys,U,x0,HessA);
      %                  Below, B is U, C is x0, and D is HessA! 
      [ Y, x ] = ldsimt( As, Bs, Cs, Ds, B, C, D );
   else
      error( 'Wrong number of input arguments' )
   end
   %
else
   % The system matrices are directly specified.
   % General call     [Y,x] = DSIMT(A,B,C,D,U,x0,HessA);
   % Special calls    [Y,x] = DSIMT(A,B,C,U,x0,HessA); 
   %                  [Y,x] = DSIMT(A,B,C,U,HessA); 
   % 
   if ni < 4,
      error( 'DSIMT needs at least 4 input parameters' )
   end
   [ m2, n2 ] = size( B );  [ m3, n3 ] = size( C );   m = n2;   l = m3;   n = n3;
   [ m4, n4 ] = size( D );
   if ni >= 5,
      m5 = size( U, 1 );
      if ni >= 7,
         [ Y, x ] = ldsimt( A, B, C, D, U, x0, HessA );
      elseif ni >= 6,
         % Special calls    [Y,x] = DSIMT(A,B,C,D,U,x0);
         %                  [Y,x] = DSIMT(A,B,C,U,x0,HessA);
         if m4 == l && n4 == m && m5 == m && length( x0 ) == n,             
            [ Y, x ] = ldsimt( A, B, C, D, U, x0 );
         else             
            % Below, D means U, and U means x0, x0 means HessA !
            if length( U ) < n,
               error( 'The initial state vector is too short' )
            end
            [ Y, x ] = ldsimt( A, B, C, zeros( l, m ), D, U, x0  );
         end
      else
         % Special calls    [Y,x] = DSIMT(A,B,C,D,U);
         %                  [Y,x] = DSIMT(A,B,C,U,x0);
         if m4 == l && n4 == m && m5 == m,             
            [ Y, x ] = ldsimt( A, B, C, D, U );
         else             
            % Below, D means U, and U means x0 !
            if length( U ) < n,
               error( 'The initial state vector is too short' )
            end
            [ Y, x ] = ldsimt( A, B, C, zeros( l, m ), D, U );
         end
      end
   else
      % Special call     [Y,x] = DSIMT(A,B,C,U);  
      if m4 == m,
         % Below, D means U !
         [ Y, x ] = ldsimt( A, B, C, zeros( l, m ), D );
      else
         error( 'The fourth input argument is wrong' )
      end
   end
end
%
% end dsimt
