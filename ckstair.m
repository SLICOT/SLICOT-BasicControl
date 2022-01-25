function result = ckstair( sys, kstair, form, triangs )
%CKSTAIR Checks that a system SYS = (A,B,C,D) is in a staircase form.
%
%        CKSTAIR(SYS,KSTAIR)  returns 1 if the system SYS is in a
%        controllability staircase form with stair sizes defined 
%        in KSTAIR, and 0 otherwise. 
%        
%        CKSTAIR(SYS,KSTAIR,FORM,TRIANGS)  returns 1 if the system SYS
%        is in a staircase form specified by KSTAIR, FORM, and TRIANGS,
%        and 0 otherwise. 
%        
%        FORM is an integer scalar specifying the staircase form to
%        be checked out, as follows:
%        FORM = 1 :  controllability staircase form;
%        FORM = 2 :  observability staircase form.
%        Default:  FORM = 1.
%
%        TRIANGS is an integer scalar specifying the form of the
%        stairs, as follows:
%        TRIANGS = 0 :  general stairs;
%        TRIANGS = 1 :  upper triangular stairs (for FORM = 1);
%        TRIANGS = 2 :  lower triangular stairs (for FORM = 2).
%        Default:  TRIANGS = 0.
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Sep-16-2003.
%        Revised Mar-02-2009.
%

%
ni = nargin;
if ni < 2,  
    error( ['Usage: CKSTAIR(SYS,KSTAIR)', sprintf('\n'), ...
            '       CKSTAIR(SYS,KSTAIR,FORM,TRIANGS)' ] );
end
%
if ni == 2,
    form    = 1;
    triangs = 0;
elseif ni == 3,
    triangs = 0;
end
%
[A,B,C,D] = ssdata( sys );
n = size( A,1 );  [p,m] = size( D );  ns = length( kstair );
if ns > n,  error( 'KSTAIR has too many elements' );  end
%
ncont = sum( kstair );
if ncont > n,  error( [ 'sum(KSTAIR) must not exceed ', int2str( n ) ] );  end
%
if n == 0,  result = 1;  return;  end
%
if ns > 0,
    irow = kstair(1);
    if ( form == 1 && irow > m ) || ( form ~= 1 && irow > p ),  
        warning( 'SLICOT:ckstair', 'KSTAIR(1) is too large' );
    end
else
    irow = 0;
end
%
if form == 1,
    if ns > 1,  irow1 = irow + kstair(2);  else  irow1 = irow;  end
    if triangs == 0,
        result = norm( B(irow+1:n,:), 1 ) == 0;  if ~result,  return;  end
        for k = 2 : ns - 1,
            icol = irow;  irow = irow1;  irow1 = irow + kstair(k+1);
            result = norm( A(irow+1:irow1,1:icol), 1 ) == 0;  
            if ~result,  return;  end
        end
    elseif triangs == 1,
        result = ( norm( B(:,1:m-irow), 1 ) + ...
                   norm( tril( B(:,m-irow+1:m), -1 ), 1 ) ) == 0;  
        if ~result,  return;  end
        icol0 = 0;
        for k = 2 : ns - 1,
            icol = irow;  irow0 = irow1;  irow1 = irow0 + kstair(k+1);
            incr = kstair(k-1) - kstair(k);  rows = irow+1 : ncont;
            result = ( norm( A(rows,icol0+1:icol0+incr), 1 ) + ...
                       norm( tril( A(rows,icol0+incr+1:icol), -1 ), 1 ) ) == 0;  
            if ~result,  return;  end
            irow = irow0;  icol0 = icol;
        end
    end
    result = norm( A(irow1+1:n,1:irow1), 1 ) == 0;  
    if ~result,  return;  end
else
    irow = 0;  if ns > 0,  irow1 = kstair(1);  else  irow1 = 0;  end;  
    icol = irow1;
    if ns > 1,  icol1 = icol + kstair(2);  else  icol1 = icol;  end
    if triangs == 0,
        result = norm( C(:,icol+1:n), 1 ) == 0;  if ~result,  return;  end
        for k = 2 : ns - 1,
            icol = icol1;  icol1 = icol + kstair(k+1);
            result = norm( A(irow+1:irow1,icol+1:n), 1 ) == 0;  
            if ~result,  return;  end
            irow = irow1;  irow1 = icol;
        end
    elseif triangs == 1,
        result = ( norm( C(1:p-irow1,:), 1 ) + ...
                   norm( triu( C(p-irow1+1:p,:), 1 ), 1 ) ) == 0;
        if ~result,  return;  end
        for k = 2 : ns - 1,
            cols = icol+1 : n;  icol = icol1;  icol1 = icol + kstair(k+1);
            incr = kstair(k-1) - kstair(k);
            result = ( norm( A(irow+1:irow+incr,cols), 1 ) + ...
                       norm( triu( A(irow+incr+1:irow1,cols), 1 ), 1 ) ) == 0;  
            if ~result,  return;  end
            irow = irow1;  irow1 = icol;
        end
    end
    result = norm( A(irow+1:ncont,ncont+1:n), 1 ) == 0;  
    if ~result,  return;  end
end
%
% end ckstair
