function [o1,o2,o3] = haconv(i1,i2,i3)
% HACONV  Converts a Hamiltonian matrix,
%
%           [ A   G ]
%      H =  [       ],
%           [ Q  -A']
%
%    where G and Q are symmetric / Hermitian matrices, between
%    various storage representations.
%
%    [A,G,Q] = HACONV(H)  extracts the n-by-n matrices A, G and
%    Q from the 2n-by-2n matrix H.
%
%    [A,QG] = HACONV(H)  stores Q and G in an n-by-n+1 array QG
%    in compressed format.
%
%    [A,G,Q] = HACONV(A,QG)
%    S = haconv(A,QG)
%    [A,QG] = HACONV(A,G,Q)
%    S = HACONV(A,G,Q)
%
%    See also SHCONV.

%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%    Based on HAPACK package, http://www.tu-chemnitz.de/mathematik/hapack/
%
%    Revisions:
%    V. Sima, June 2008, March 2009.
%

ni = nargin;  no = nargout;
%
if ni < 1,
   error( 'Not enough input arguments' );
end
if no > 3,
   error( 'Too many output arguments' );
end
%
if ni == 1,
   n = floor( size( i1, 1 ) / 2 );
   if ( size( i1, 1 ) ~= 2*n ) || ( size( i1, 2 ) ~= 2*n ),
      error( 'H is not a 2n-by-2n matrix' );
   end
   if no <= 1,
      o1 = i1;
   elseif no == 2,
      o1 = i1(1:n,1:n);
      o2 = zeros( n, n+1 );
      o2(:,1:n)   = tril( i1(n+1:2*n,1:n) );
      o2(:,2:n+1) = o2(:,2:n+1) + triu( i1(1:n,n+1:2*n) );
   elseif no == 3,
      o1 = i1(1:n,1:n);
      o2 = i1(1:n,n+1:2*n);
      o3 = i1(n+1:2*n,1:n);
   end
elseif ni == 2,
   n = size( i1, 1 );
   if ( size( i1, 2 ) ~= n ),
      error( 'A is not a square matrix' );
   end
   if ( size( i2, 1 ) ~= n ) || ( size( i2, 2 ) ~= n+1 ),
      error( 'QG is not an n-by-(n+1) matrix' );
   end
   if no <= 1,
      o1 = [ i1, triu( i2(:,2:n+1) )+triu( i2(:,2:n+1),1 )'; ...
             tril( i2(:,1:n) )+tril( i2(:,1:n),-1 )', -i1' ];
   elseif no == 2,
      o1 = i1;
      o2 = i2;
   elseif no == 3,
      o1 = i1;
      o2 = triu( i2(:,2:n+1) ) + triu( i2(:,2:n+1),1 )';
      o3 = tril( i2(:,1:n) )   + tril( i2(:,1:n),-1 )';
   end
elseif ni == 3,
   n = size( i1, 1 );
   if ( size( i1, 2 ) ~= n ),
      error( 'A is not a square matrix' );
   end
   if ( size( i2, 1 ) ~= n ) || ( size( i2, 2 ) ~= n ),
      error( 'G is not an n-by-n matrix' );
   end
   if ( size( i3, 1 ) ~= n ) || ( size( i3, 2 ) ~= n ),
      error( 'Q is not an n-by-n matrix' );
   end
   if no <= 1,
      o1 = [i1, i2; i3, -i1'];
   elseif no == 2,
      o1 = i1;
      o2 = zeros( n, n+1 );
      o2(:,1:n)   = tril( i3 );
      o2(:,2:n+1) = o2(:,2:n+1) + triu( i2 );
   elseif no == 3,
      o1 = i1;
      o2 = i2;
      o3 = i3;
   end
else
   error( 'Too many input arguments.' );
end
