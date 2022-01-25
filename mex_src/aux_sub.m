function aux_sub
%
% Function for compiling auxiliary subroutines for building MEX-executables using
% the MATLAB built-in FORTRAN compiler (Linux version, 64 bit).

% Contributor:
% V. Sima, Jan. 2022.
%
% Revisions:
% -

flags = 'FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer -fdefault-integer-8" -largeArrayDims';
%
%% Set location of the source files for auxiliary subroutines.
%
aux_src = '';
%
%% Compiling.
%
aux_sub = {
    'MA02KD', ...
    'MA02KI', ...
    'MA02KS', ...
    'MA02KV', ...
    'MA02KW', ...
    'MA02KZ', ...
    'MA02LD', ...
    'MA02LZ', ...
    'MA02QD', ...
    };
obj = [];
for k = 1:length(aux_sub)
    file = aux_sub{k};
    fprintf( 'mex -c %s %s%s.f\n', flags, aux_src, file );
    eval( sprintf( 'mex -c %s %s%s.f', flags, aux_src, file ) );
    obj = [obj ' ' file '.o'];
end
%
% Archive the object files.
%
unix('ar r aux_sub.a *.o');
unix('rm *.o');
