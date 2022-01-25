function [Us,Uc,w] = sHHstab(A,DE,B,FG,meth)
% SHHSTAB  Complete stable right deflating subspace of a skew-Hamiltonian/
%          Hamiltonian matrix pencil. The stable companion subspace can
%          also be returned for a pencil with the skew-Hamiltonian matrix 
%          in factored form.
%
%    US = SHHSTAB(S,H) is an orthonormal basis of the stable right deflating
%    subspace of a 2*n-by-2*n skew-Hamiltonian/Hamiltonian matrix pencil
%    lambda*S - H, where S is skew-Hamiltonian and H is Hamiltonian,
%
%          (  A  D  )         (  B  F  )
%      S = (        ),    H = (        ),
%          (  E  A' )         (  G -B' )
%
%    with D and E skew-symmetric/skew-Hermitian matrices, and F and G 
%    symmetric/Hermitian matrices. 
%
%    US = SHHSTAB(A,DE,B,FG) requires A and DE, and B and FG, contain the
%    matrices S and H, respectively, in compressed format. 
%    If D and E are available instead, DE can be obtained using
%      [~,DE] = shconv(A,D,E);
%    Similarly, if F and G are available, FG can be obtained using
%      [~,FG] = haconv(B,F,G);
%
%    US = SHHSTAB(S,B,FG), or
%    US = SHHSTAB(A,DE,H) require that either H or S is in compressed format.
%
%    US = SHHSTAB(Z,H,'factor'), or
%    US = SHHSTAB(Z,B,FG,'factor') require that S is factored, S = J'*Z'*J*Z
%    (with J = [0 I; -I 0]), and the factor Z is given.
%
%    [US,UC] = SHHSTAB(Z,H,'factor'), or
%    [US,UC] = SHHSTAB(Z,B,FG,'factor') also return an orthonormal basis
%    of the stable companion subspace.
%
%    [US,w]    = SHHSTAB(...), or 
%    [US,UC,w] = SHHSTAB(...), also return the eigenvalues.
%
%    [...]    = SHHSTAB(...,'SVD') require that singular value decomposition
%    be used for finding the subspace(s). By default, QR decomposition with
%    column pivoting is used.
%
%    Note that SHHSTAB supports both real and complex skew-Hamiltonian/
%    Hamiltonian matrix pencils.
%
%    See also EIG, SHHEIG.

%    This is a MATLAB gateway routine to the SLICOT routines
%    MB03FZ, MB03LF, MB03LD, and MB03LZ.
%
%    RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%    Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%    Contributor
%    V. Sima, Nov. 2012.
%
%    Revision
%    -

if nargin < 2,
   error('At least two input arguments required.');
end
%
errMsg1 = 'String argument is an unknown option';
errMsg2 = 'argument should be string';
errMsg3 = 'is not an n-by-n matrix';
errMsg4 = 'is not a 2n-by-2n matrix';
errMsg5 = 'should be in a compressed format';
errMsg6 = 'matrices should have the same row size';
%
is_fact = false;
is_SVD = false;
%
isr = isreal(A) && isreal(DE);
if nargin == 2,
   compressed = false;
elseif nargin == 3,
   if ischar(B),
      is_fact = isequal(lower(B(1)),'f');
      is_SVD  = isequal(lower(B(1)),'s');
      if ~is_fact && ~is_SVD,  error(errMsg1);  end
      compressed = false;
   else
      compressed = true;
      isr = isr && isreal(B);
   end
elseif nargin == 4,
   if ischar(B),
      is_fact = isequal(lower(B(1)),'f');
      if ~is_fact,  error(errMsg1);  end
      if ischar(FG),
         is_SVD = isequal(lower(FG(1)),'s');
         if ~is_SVD,  error(errMsg1);  end
      else
         error(['The fourth ' errMsg2]);
      end
      compressed = false;
   elseif ischar(FG),
      is_fact = isequal(lower(FG(1)),'f');
      is_SVD  = isequal(lower(FG(1)),'s');
      if ~is_fact && ~is_SVD,  error(errMsg1);  end
      compressed = true;
      isr = isr && isreal(B);
   else
      compressed = true;
      isr = isr && isreal(B) && isreal(FG);
   end
elseif nargin == 5,
   isr = isr && isreal(B);
   if ischar(FG),
      is_fact = isequal(lower(FG(1)),'f');
      if ~is_fact,  error(errMsg1);  end
   else
      isr = isr && isreal(FG);
   end
   if ~ischar(meth),  error(['The fifth ' errMsg2]);  end
   is_SVD = isequal(lower(meth(1)),'s');
   if ~is_SVD,  error(errMsg1);  end
   compressed = true;
else
   error('Too many input arguments');
end
%
if is_fact 
   mxout = 3;
   if nargout > 1,  compu = 1;  else  compu = 0;  end
else
   mxout = 2;
end
%
if nargout > mxout,
   error('Too many output arguments');
end
%
if ~compressed,
   n = fix(size(A,1)/2);
   if (size(A,1) ~= 2*n) || (size(A,2) ~= 2*n),
      if is_fact,
         error(['Z ' errMsg4]);
      else
         error(['S ' errMsg4]);
      end
   elseif (size(DE,1) ~= 2*n) || (size(DE,2) ~= 2*n),
      error(['H ' errMsg4]);
   end
   [B,FG] = haconv(DE);
   if ~is_fact,  [A,DE] = shconv(A);  end
elseif nargin == 3 || ( nargin == 4 && is_SVD ),
   [md,ld] = size(DE);  [mf,lf] = size(B);
   if ld == md+1,
      compressedDE = true;  n = md;
      compressedFG = false;
      if lf == mf+1,  error('Wrong dimensions');  end 
   elseif lf == mf+1,
      compressedDE = false;
      compressedFG = true;  n = mf;
   else
      error(['Either S or H ' errMsg5]);
   end
elseif nargin == 4,
   if is_fact,
      [mf,lf] = size(B);
      if lf == mf+1,
         compressedDE = false;
         compressedFG = true;  n = mf;
      end
   else
      [md,ld] = size(DE);  [mf,lf] = size(FG);
      if ld ~= md+1 || lf ~= mf+1,  error(['Both S and H ' errMsg5]);  end
      compressedDE = true;  n = md;
      compressedFG = true;
      if md ~= n || mf ~= n,  error(['DE and FG ' errMsg6]);  end
   end
else
   if is_fact,
      n = size(A,1)/2;
      compressedDE = false;
      compressedFG = true;
      [mf,lf] = size(B);
      if lf ~= mf+1,  error(['H ' errMsg5]);  end
   else
      n = size(A,1);
      compressedDE = true;
      compressedFG = true;
      [md,ld] = size(DE);  [mf,lf] = size(FG);
      if ld ~= md+1 || lf ~= mf+1,  error(['Both S and H ' errMsg5]);  end
      if md ~= n || mf ~= n,  error(['All ' errMsg6]);  end
   end
end
%
if compressed,
   if compressedDE,
      if (size(A,1) ~= n),
         error(['A and DE ' errMsg6]);
      elseif (size(A,2) ~= n),
         error(['A ' errMsg3]);
      elseif ~compressedFG,
         if (size(B,1) ~= 2*n) || (size(B,2) ~= 2*n),  error(['H ' errMsg4]);  end
         [B,FG] = haconv(B);
      elseif (size(B,1) ~= n) || (size(B,2) ~= n),
         error(['B ' errMsg3]);
      elseif (size(FG,1) ~= n) || (size(FG,2) ~= n+1),
         error('FG is not an n-by-(n+1) matrix');
      end
   elseif compressedFG,
      if (size(A,1) ~= 2*n) || (size(A,2) ~= 2*n),
         if is_fact,
            error(['Z ' errMsg3]);
         else
            error(['S ' errMsg3]);
         end
      elseif ~is_fact && ( (size(DE,1) ~= n) || (size(DE,2) ~= n) ),
         error(['B ' errMsg3]);
      end
      if is_fact,  
         if ischar(B),  [B,FG] = haconv(DE);  else  FG = B;  B = DE;  end
      else
         FG = B;  B = DE;
         [A,DE] = shconv(A);
      end
   end
end
%
compq = 1;
if is_SVD,  orthm = 2;  else  orthm = 1;  end
%
if isr,
   if is_fact,
      [wr,wi,bet,Us,Uc] = skewHamildeflf(A,B,FG,compq,compu,orthm);
   else
      [wr,wi,bet,Us] = skewHamildefl(A,DE,B,FG,compq,orthm);
   end
   if nargout == mxout,
      if numel(bet) > n,
         w = (wr+1i*wi)./bet(1:n).*bet(2*n+1).^bet(n+1:2*n);
      else
         w = (wr+1i*wi)./bet;
      end
      w = [w;-w];
   end
else
   if is_fact,
      [x,bet,Us,Uc] = skewHamildeflfZ(A,B,FG,compq,compu,orthm);
   else
      [x,bet,Us] = skewHamildeflZ(A,DE,B,FG,compq,orthm);
   end
   if nargout == mxout,  w = x./bet;  end
end
if nargout == mxout,  
   [~,ii] = sort(real(w));
   w = w(ii);
   if nargout == 2,  Uc = w;  end
end

