function[u,beta]=house(x)
% [u, beta] = house(x)
%
% input:
%
%    x - real n-vector
%
% output:
%   
%    u - real n-vector containing the Householder vector
% beta - the real scalar factor of the Householder transformation
%
% purpose:
%
%    To compute a Householder transformation which transforms the
%    real n-vector x into the first unit vector.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Revisions: March 2009

m=max(abs(x));
if m==0,
  beta=0;u=x;
else
  x=x/m;u=x;

  hh=sign(u(1));
  if u(1)==0,
    hh=1;
  end
  uu= norm(u,2);

  sigma=hh*uu;

  u(1)=u(1)+sigma;
  beta=2/(u'*u); 

end
