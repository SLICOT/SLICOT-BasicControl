function w = sHHeig(A,DE,B,FG,flag)
% SHHEIG  Eigenvalues of a skew-Hamiltonian/(skew-)Hamiltonian matrix pencil.
%
%    w = SHHEIG(S,H) is a vector containing the eigenvalues of a
%    2*n-by-2*n skew-Hamiltonian/Hamiltonian matrix pencil lambda*S - H,
%    where S is skew-Hamiltonian and H is Hamiltonian,
%
%          (  A  D  )         (  B  F  )
%      S = (        ),    H = (        ),
%          (  E  A' )         (  G -B' )
%
%    with D and E skew-symmetric/skew-Hermitian matrices, and F and G 
%    symmetric/Hermitian matrices. 
%
%    w = SHHEIG(A,DE,B,FG) requires A and DE, and B and FG, contain the
%    matrices S and H, respectively, in compressed format. 
%    If D and E are available instead, DE can be obtained using
%      [~,DE] = shconv(A,D,E);
%    Similarly, if F and G are available, FG can be obtained using
%      [~,FG] = haconv(B,F,G);
%
%    w = SHHEIG(S,B,FG), or
%    w = SHHEIG(A,DE,H) require that either H or S is in compressed format.
%
%    w = SHHEIG(Z,H,'factor'), or
%    w = SHHEIG(Z,B,FG,'factor') require that S is factored, S = J'*Z'*J*Z
%    (with J = [0 I; -I 0]), and the factor Z is given.
%
%    w = SHHEIG(S,H,'skew2'),    or 
%    w = SHHEIG(S,B,FG,'skew2'), or
%    w = SHHEIG(A,DE,H,'skew2'), or 
%    w = SHHEIG(A,DE,B,FG,'skew2'), compute the eigenvalues of a real
%    skew-Hamiltonian/skew-Hamiltonian matrix pencil, with
%
%          (  A  D  )         (  B  F  )
%      S = (        ),    H = (        ),
%          (  E  A' )         (  G  B' )
%
%    where D, E, F, and G are skew-symmetric matrices. 
%
%    w = SHHEIG(Z,H,'factor','skew2'), or 
%    w = SHHEIG(Z,B,FG,'factor','skew2') require that S is factored, and
%    its factor Z is given.
%
%    Note that SHHEIG supports both real and complex skew-Hamiltonian/
%    Hamiltonian matrix pencils, but only real skew-Hamiltonian/
%    skew-Hamiltonian matrix pencils.
%
%    See also EIG, SHHSTAB

%    This is a MATLAB gateway routine to the SLICOT routines
%    MB04AD, MB04AZ, MB04BD, MB04BZ, MB04ED, and MB04FD.
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
if ( nargout > 1 ),
   error('Too many output arguments');
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
is_sHsH = false;
%
isr = isreal(A) && isreal(DE);
if nargin == 2,
   compressed = false;
elseif nargin == 3,
   if ischar(B),
      is_fact = isequal(lower(B(1)),'f');
      is_sHsH = isequal(lower(B(1)),'s');
      if ~is_fact && ~is_sHsH,  error(errMsg1);  end
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
         is_sHsH = isequal(lower(FG(1)),'s');
         if ~is_sHsH,  error(errMsg1);  end
      else
         error(['The fourth ' errMsg2]);
      end
      compressed = false;
   elseif ischar(FG),
      is_fact = isequal(lower(FG(1)),'f');
      is_sHsH = isequal(lower(FG(1)),'s');
      if ~is_fact && ~is_sHsH,  error(errMsg1);  end
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
   if ~ischar(flag),  error(['The fifth ' errMsg2]);  end
   is_sHsH = isequal(lower(flag(1)),'s');
   if ~is_sHsH,  error(errMsg1);  end
   compressed = true;   
else
   error('Too many input arguments');
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
   if is_sHsH,
      [B,FG] = shconv(DE);
   else
      if ~isr || ~is_fact,  [B,FG] = haconv(DE);  end
   end
   if ~is_fact,  [A,DE] = shconv(A);  end
elseif nargin == 3 || ( nargin == 4 && ( is_fact || is_sHsH ) ),
   [md,ld] = size(DE);  [mf,lf] = size(B);
   if ld == md+1,
      compressedDE = true;  n = md;
      if lf == mf+1,  error('Wrong dimensions');  end 
      compressedFG = false;
   elseif lf == mf+1,
      compressedDE = false;
      compressedFG = true;  n = mf;
   else
      error(['Either S or H ' errMsg5]);
   end
elseif nargin == 4,
   [md,ld] = size(DE);  [mf,lf] = size(FG);
   if ld ~= md+1 || lf ~= mf+1,  error(['Both S and H ' errMsg5]);  end
   compressedDE = true;  n = md;
   compressedFG = true;
   if md ~= n || mf ~= n,  error(['DE and FG ' errMsg6]);  end
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
         if is_sHsH,
            [B,FG] = shconv(B);
         else
            [B,FG] = haconv(B);
         end
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
         if isr && ~is_sHsH,  DE = haconv(DE,B);  end
         if ischar(B),  [B,FG] = shconv(DE);  else  FG = B;  B = DE;  end
      else
         FG = B;  B = DE;
         [A,DE] = shconv(A);
      end
   end
end
%
if isr,
   if is_sHsH,
      if is_fact,
         [wr,wi,bet] = skewHamil2feig(A,B,FG);
      else 
         [wr,wi,bet] = skewHamil2eig(A,DE,B,FG);
      end
      w = (wr+1i*wi)./bet;
      w = [w;w];
   else
      if is_fact,
         [wr,wi,bet] = symplURV(A,DE);
      else 
         [wr,wi,bet] = skewHamileig(A,DE,B,FG);
      end
      if numel(bet) > n,
         w = (wr+1i*wi)./bet(1:n).*bet(2*n+1).^bet(n+1:2*n);
      else
         w = (wr+1i*wi)./bet;
      end
      w = [w;-w];
   end
else
   if is_sHsH
      error(['Problems with complex skew-Hamiltonian/skew-Hamiltonian ', ...
             'matrix pencils cannot be solved.'])
   elseif is_fact,
      [x,bet] = symplURVZ(A,B,FG);
   else
      [x,bet] = skewHamileigZ(A,DE,B,FG);
   end
   w = x./bet;
end
[~,ii] = sort(real(w));
w = w(ii);

