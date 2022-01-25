function p = slpole( sys, E )
%SLPOLE  Computes the poles of a standard or descriptor system
%        SYS = (A,B,C,D,E), with square matrices A and E.
%
%        P = SLPOLE(SYS)  computes the poles of the system SYS. 
%        
%        P = SLPOLE(A)    computes the eigenvalues of the  
%        matrix A.
%        
%        P = SLPOLE(A,E)  computes the generalized eigenvalues of 
%        the matrix pencil (A,E).
%        
%        Note: This function cannot work with a singular matrix E,
%              or with rectangular A and E. Use POLEZERO(Z) instead.
%
%        See also POLEZERO, POLEZEROZ
%

%        RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Dec-22-2003.
%        Revised Jul-14-2004; Dec-13-2008, Mar-03-2009.
%

%
ni = nargin;
if ni < 1,  
    error( ['Usage: P = SLPOLE(SYS)', sprintf('\n'), ...
            '       P = SLPOLE(A)',   sprintf('\n'), ...
            '       P = SLPOLE(A,E)' ] );
end
%
if ni == 1,
    if isa( sys, 'ss' ),
        E = sys.e;
        iscmplx = any( any( imag( sys.a ) ) ) | any( any( imag( E ) ) );
        if iscmplx,
            if ~isempty( E ) && ~isequal( E, eye( size( sys.a ) ) ),
                [ poles, betap ] = polezeroz( 0, sys.a, sys.e );
                p = poles ./ betap;
            else  
                p = polezeroz( 0, sys.a );
            end
        else
            if ~isempty( E ) && ~isequal( E, eye( size( sys.a ) ) ),
                [ rpoles, ipoles, betap ] = polezero( 0, sys.a, sys.e );
                p = complex( rpoles, ipoles ) ./ betap;
            else
                [ rpoles, ipoles ] = polezero( 0, sys.a );
                p = complex( rpoles, ipoles );
            end
        end
    else
        iscmplx = any( any( imag( sys ) ) );
        if iscmplx,
            p = polezeroz( 0, sys );
        else
            [ rpoles, ipoles ] = polezero( 0, sys );
            p = complex( rpoles, ipoles );
        end
    end
elseif ni == 2,
    iscmplx = any( any( imag( sys ) ) ) | any( any( imag( E ) ) );
    if iscmplx,
        if ~isequal( E, eye( size( sys ) ) ),  
            [ poles, betap ] = polezeroz( 0, sys, E );
            p = poles ./ betap;
        else
            p = polezeroz( 0, sys );
        end
    else
        if ~isequal( E, eye( size( sys ) ) ),
            [ rpoles, ipoles, betap ] = polezero( 0, sys, E );
            p = complex( rpoles, ipoles ) ./ betap;
        else
            [ rpoles, ipoles ] = polezero( 0, sys );
            p = complex( rpoles, ipoles );
        end
    end
end
%
% end slpole
