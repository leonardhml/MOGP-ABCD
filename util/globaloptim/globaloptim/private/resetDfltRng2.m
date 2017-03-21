function resetDfltRng2(rngstate)
%RESETDFLTRNG2 Reset the random number generator to a saved state.
%
%   This function is private to the GADS solvers

%   Copyright 2008-2010 The MathWorks, Inc.

if ~isempty(rngstate)
    if isfield(rngstate,'state') && isfield(rngstate,'type')
        dflt = RandStream.getGlobalStream;
        if isempty(rngstate.type) || isequal(rngstate.type,dflt.Type)
            try
                dflt.State = rngstate.state;
            catch ME
                error(message('globaloptim:resetDfltRng2:InvalidRNGState'));
            end
        else
            error(message('globaloptim:resetDfltRng2:RNGType', rngstate.type, dflt.Type));
        end
    else
        error(message('globaloptim:resetDfltRng2:InvalidRNGField'));
    end
end
