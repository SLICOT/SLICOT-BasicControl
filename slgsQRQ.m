function [sysQR,U] = slgsQRQ(sys,job,compu,U1)
%SLGSQRQ Transform the pair (A,E) of a descriptor system to a
%        QR- or RQ-coordinate form.
%
%        SYSQR = SLGSQRQ(SYS)  applies an orthogonal transformation
%        matrix U to a descriptor system pair (A-lambda E,B) such that
%        the transformed descriptor system pair (in SYSQR),
%        (U'*A-lambda U'*E, U'*B), has the descriptor matrix U'*E in an
%        upper trapezoidal form. This is the QR-coordinate form.
%        The matrix U is not computed.
%
%        [SYSQR,U] = SLGSQRQ(SYS,JOB,COMPU,U1)  has additional input
%        arguments.
%
%        JOB is an integer scalar specifying whether the QR- or 
%        RQ-coordinate form is desired, as follows:
%        JOB = 0 :  QR-coordinate form above;
%        JOB = 1 :  RQ-coordinate form, i.e., an orthogonal transformation
%        matrix U is applied to the descriptor system pair (C,A-lambda E)
%        such that the transformed descriptor system pair (in SYSQR),
%        (C*U,A*U-lambda E*U), has the descriptor matrix E*U in an upper
%        trapezoidal form.
%        Default:  JOB = 0.
%
%        COMPU is an integer scalar indicating what should be done
%        with matrix U, as follows:
%        COMPU = 0 :  do not compute U;
%        COMPU = 1 :  U is initialized to the unit matrix, and
%                     the orthogonal matrix U is returned;
%        COMPU = 2 :  U is initialized to an orthogonal matrix U1
%                     and the product U1*U is returned.
%        Default:  COMPU = 0.
%        If COMPU = 0, U must not be an output argument.
%
%        Note: This function cannot work with a singular matrix E. 
%              Use GSYSTRA instead.
%        Note that GSYSTRA can work with rectangular matrices A and E.
%
%        See also SLGSBAL, SLGSHES, SLGSRSF, SLGSSVD, GSYSTRA
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 30-04-2003.
%
%        Revisions: 03-03-2009.

ni = nargin;
no = nargout;
%
[A,B,C,D,E] = dssdata( sys );
if norm( E - eye( size( E ) ), 1 ) == 0,
   % Standard system: nothing to do.
   sysQR = sys;
   if no == 2,  U = eye( size( A ) );  end
%
else
   % Descriptor system.
   if ni <= 1,  job   = 0;  end
   if ni <= 2,  compu = 0;  end
   %
   if job   < 0 || job   > 1,  error('Wrong value for JOB');    end
   if compu < 0 || compu > 2,  error('Wrong value for COMPU');  end
   %
   task = job + 2;
   if job == 0,
      if compu == 0, 
         [A,E,B] = gsystra(task,A,E,B);
      elseif compu == 1,
         [A,E,B,U] = gsystra(task,A,E,B,compu);
      else
         [A,E,B,U] = gsystra(task,A,E,B,compu,U1);
      end
   else
      if compu == 0, 
         [A,E,C] = gsystra(task,A,E,C);
      elseif compu == 1,
         [A,E,C,U] = gsystra(task,A,E,C,compu);
      else
         [A,E,C,U] = gsystra(task,A,E,C,compu,U1);
      end
   end
   sysQR = dss( A,B,C,D,E,sys );
end
%
% end slgsQRQ
