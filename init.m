% Load the random number generator seed.

%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   V. Sima, Katholieke Univ. Leuven, Belgium, June 1999.
%
%   Revisions:
%   V. Sima, May 2005, March 2009.
%   
load Seed
vrs = version;
if str2double( vrs(1) ) >= 7 && str2double( vrs(3) ) >= 7,
    defaultStream = RandStream.getDefaultStream;  
    defaultStream.State = savedState;
else
    rand( 'state', seed );
end
