function [state, optnew,optchanged] = gaoutputoptions2(options,state,flag)
%GAOUTPUTOPTIONS2 Display non default options settings.
%   STATE = GAOUTPUTOPTIONS2(OPTIONS,STATE,FLAG) is used as a output
%   function for displaying which options parameters are not default.
%
%   Example:
%    Create an options structure that uses GAOUTPUTOPTIONS2
%    as the output function
%     options = optimoptions2('ga2','OutputFcn',@gaoutputoptions2)
% 
%    (Note: If calling gamultiobj2, replace 'ga2' with 'gamultiobj2') 

%   Copyright 2003-2015 The MathWorks, Inc.
 
% use the code generation tool to get the text to display
optnew = options;
optchanged = false;
if(strcmp(flag,'init'))
    code = generateCode2(options)
end