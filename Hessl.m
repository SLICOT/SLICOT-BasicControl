function X = Hessl( H, B, tran )
%HESSL   Solves a set of systems of linear equations with an
%        upper Hessenberg coefficient matrix H.
% 
%        X = HESSL(H,B)  solves the set of systems of linear equations
%        H*X = B, where H is an upper Hessenberg matrix of order N
%        and B is an N-by-P matrix.
%
%        X = HESSL(H,B,TRAN)  has an additional input argument:
%
%        TRAN specifies whether the matrix H or its (conjugate) transpose
%        should be used.
%        TRAN = 0:  use the matrix H;
%        TRAN = 1:  use the matrix H';
%        TRAN = 2:  use the matrix H**H.
%        Default:  TRAN = 0.
%
%        See also HESSOL, HESSCOND
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Dec. 2002.
%        Revised: Oct. 2004, Mar. 2009.
%

ni = nargin;
%
if ni < 2,
    error('Usage: X = HESSL(H,B,TRAN)')
elseif ni == 2,
    tran = 0;
end
%
rH = isreal( H );  rB = isreal( B );  normh = 1;
if rH == 1 && rB == 1,
    mtype = 0;
    X = Hessol(1,H,mtype,normh,B,tran);
elseif rH == 1,
    mtype = 1;  Hc = H + 1i*zeros( size( H ) );
    X = Hessol(1,Hc,mtype,normh,B,tran);
elseif rB == 1,
    mtype = 1;  Bc = B + 1i*zeros( size( B ) );
    X = Hessol(1,H,mtype,normh,Bc,tran);
else
    mtype = 1;
    X = Hessol(1,H,mtype,normh,B,tran);
end
%
% end Hessl
