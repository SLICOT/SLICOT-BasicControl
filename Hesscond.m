function rcnd = Hesscond( H, normh )
%HESSCOND Computes an estimate of the reciprocal of the condition
%        number of an upper Hessenberg matrix H.
% 
%        RCND = HESSCOND(H)  computes an estimate of the reciprocal
%        of the condition number of H in the 1-norm.
%
%        RCND = HESSCOND(H,NORMH)  has an additional input argument:
%
%        NORMH specifies whether the 1-norm or the infinity-norm 
%        reciprocal condition number is required.
%        NORMH = 1:  1-norm;
%        NORMH = 2:  infinity-norm.
%        Default:  NORMH = 1.
%
%        See also HESSOL, HESSL
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Dec. 2002.
%        Revised: Oct. 2004, Mar. 2009.
%

ni = nargin;
%
if ni < 1,
    error('Usage: RCND = HESSCOND(H,NORMH)')
elseif ni == 1,
    normh = 1;
end
%
if isreal( H ) == 1,
    mtype = 0;
else
    mtype = 1;
end
rcnd = Hessol(2,H,mtype,normh);
%
% end Hesscond
