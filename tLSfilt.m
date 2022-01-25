function [Ef,Xf,Eb,Cs,Sn,Err,Salph,YQ] = tLSfilt( Si,Ef,Xf,Eb,Cs,Sn,Lambda,YQ )
%TLSFILT Computes a combined measurement and time update of
%        the recursive least-squares filter. It is intended
%        for testing the gateway KFILTUPD for task = 3.
%
%        [EF,XF,EB,CS,SN,ERR,SALPH] = TLSFILT(SI,EF,XF,EB,CS,SN)
%        updates the values of the input parameters, and returns 
%        other useful results. See below.
%        
%        [EF,XF,EB,CS,SN,ERR,SALPH,YQ] = TLSFILT(SI,EF,XF, ...
%                                           EB,CS,SN,LAMBDA,YQ ) 
%        has additional input and output parameters.
%        
%        Input parameters:
%        SI is a vector with one or two elements.     
%        SI(1) = Ui contains the input sample at instant i.
%        If SI has length 2, SI(2) = Yi contains the reference sample
%        at instant i. Otherwise, Yi is not used.
%        (The situation just before and just after the call of the
%        function are denoted by instant (i-1) and instant i,
%        respectively.)
%
%        EF is the square root of exponentially weighted forward
%        prediction error energy at instant (i-1).  EF >= 0.0.
%
%        XF is an L-vector containing the transformed forward prediction
%        variables at instant (i-1).
%
%        EB is an (L+1)-vector; the leading L elements must contain the
%        normalized a posteriori backward prediction error residuals
%        of orders zero through L-1, respectively, at instant (i-1),
%        and EB(L+1) must contain the square-root of the so-called
%        "conversion factor" at instant (i-1).
%
%        CS is an L-vector containing the cosines of the rotation angles used
%        in time updates, at instant (i-1).
%
%        SN is an L-vector containing the sines of the rotation angles used
%        in time updates, at instant (i-1).
%
%        LAMBDA is an optional real scalar specifying the square root of the
%        forgetting factor. For tracking capabilities and exponentially
%        stable error propagation, LAMBDA < 1.0 (strict inequality) should
%        be used.  0 < LAMBDA <= 1.
%        Default :  LAMBDA = 1.
%
%        YQ is an optional L-vector. If length(SI) = 2, YQ contains the
%        orthogonally transformed reference vector at instant (i-1).
%        These elements are also the tap multipliers of an equivalent
%        normalized lattice least-squares filter.
%        If length(SI) = 1, this parameter is not used.
%
%        Output parameters:
%        EF is the square root of the exponentially weighted forward
%        prediction error energy at instant i.
%
%        XF is L-vector containing the transformed forward prediction
%        variables at instant i.
%
%        EB is an (L+1)-vector containing the normalized a posteriori
%        backward prediction error residuals, plus the square root
%        of the conversion factor at instant i.
%
%        CS is an L-vector containing the cosines of the rotation angles
%        at instant i.
%
%        SN is an L-vector containing the sines of the rotation angles
%        at instant i.
%
%        ERR is a vector with one or two elements.     
%        ERR(1) = Ep contains the a posteriori forward prediction
%        error residual.
%        If length(SI) = 2, ERR(2) = Eo contains the a posteriori output
%        error residual from the least-squares filter at instant i.
%
%        SALPH is an L-vector: the element SALPH(i), i=1,...,L, contains
%        the opposite of the i-(th) reflection coefficient for the
%        least-squares normalized lattice predictor (whose value is
%        -SALPH(i)).
%
%        YQ is an optional L-vector. If length(SI) = 2, it contains
%        the orthogonally transformed reference vector at instant i.
%        If length(SI) = 1, YQ must not be specified as an output
%        parameter.
%
%        See also KFILTUPD
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Nov-18-2003.
%        Revised Jul-14-2004, Mar-03-2009.
%

ni = nargin;
if ni < 6 || nargout == 0,  
    error( ['Usage: [EF,XF,EB,CS,SN,ERR,SALPH]    = TLSFILT(SI,EF,XF,EB,CS,SN)', sprintf('\n'), ...
            '       [EF,XF,EB,CS,SN,ERR,SALPH,YQ] = TLSFILT(SI,EF,XF,EB,CS,SN,LAMBDA,YQ)' ] );
end
%
if ni < 7,  Lambda = 1;  end
%
% Forward prediction rotations.
%
L = numel( Xf );
fn = Si(1);  
for k = 1 : L
   Xfk   = Xf(k)*Lambda;
   Xf(k) = Sn(k)*fn + Cs(k)* Xfk;
   fn    = Cs(k)*fn - Sn(k)* Xfk;
end
Err = fn*Eb(L+1);
%
% Update the square root of the prediction energy.
%
Ef = Ef*Lambda;
tmp = norm( [ fn, Ef ] );
if ( tmp < eps ),  fn = 0;  else  fn = fn*Eb(L+1)/tmp;  end
Ef = tmp;
%
% Calculate the reflection coefficients and the backward prediction errors.
%
Salph = zeros( 1,L );
for k = L : -1 : 1
   G = givens( tmp, Xf(k) );  tmp = norm( [ tmp, Xf(k) ] );
   Eb(k+1)  = G(2,:)*[fn; Eb(k)];  fn = G(1,:)*[fn; Eb(k)];  
   Salph(k) = G(1,2);
end
Eb(1) = fn;
%
% Update to new rotation angles.
%
nrm = norm( Eb(1:L) );  tmp = sqrt( ( 1 + nrm )*( 1 - nrm ) );  Eb(L+1) = tmp;
%
for k = L : -1 : 1
   G = givens( tmp, Eb(k) );  tmp = norm( [ tmp, Eb(k) ] );
   Cs(k) = G(1,1);  Sn(k) = G(1,2);
end
%
% Joint process section.
%
if numel( Si ) == 2,
   fn = Si(2);
   for k = 1 : L
      YQk   = YQ(k)*Lambda;
      YQ(k) = Sn(k)*fn + Cs(k)*YQk;
      fn    = Cs(k)*fn - Sn(k)*YQk;
   end
   Err(2) = fn*Eb(L+1);
end
%
% end tLSfilt

