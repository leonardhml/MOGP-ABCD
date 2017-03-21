function typeValueChecker2(type,value,propertyName,valueData)
%TYPEVALUECHECKER2 Checks the value entered with respect to type.
%
%   msgstruct = TYPEVALUECHECKER2(type,value,propertyName) returns an error
%   message structure containing the fields message and identifier. An
%   appropriate error message and error ID, which includes the propertyName
%   is returned in the structure. If the value is of specified type
%   msgstruct will be an empty message structure.
%
%   msgstruct = TYPEVALUECHECKER2(type,value,propertyName,valueData) is used
%   for certain types (e.g. stringsType, boundedReal) to check the values
%   against the possible ones given in valueData. For stringsType valueData
%   is a cell array of strings. For boundedReal it is a 1x2 vector of
%   double.
%   type - one of the following: 'displayType', 'nonNegReal', 'posReal',
%   'posInteger', 'stringsType', 'boundedReal', 'functionOrCellArray'.

%   Copyright 2009-2011 The MathWorks, Inc.

if nargin < 3
    propertyName = '';
end

switch type
    case {'displayType'}
        % One of these strings: on, off, none, iter, final
        valueValid =  ischar(value) && any(strcmp(value, ...
            {'on';'off';'none';'iter';'final'}));
        if ~valueValid
            errid = 'globaloptim:typeValueChecker2:NotADisplayType';
            error(message(errid,'Display','off','on','iter','final'));
        end
    case {'nonNegReal'}
        valueValid = isa(value,'double') && isreal(value) && isscalar(value) ...
            && (value >= 0);
        if ~valueValid
            errid = 'globaloptim:typeValueChecker2:NotANonNegReal';
            error(message(errid,propertyName));
        end
    case {'posReal'}
        valueValid = isa(value,'double') && isreal(value) && isscalar(value) ...
            && (value > 0);
        if ~valueValid
            errid = 'globaloptim:typeValueChecker2:NotAPosRealNum';
            error(message(errid,propertyName));
        end        
    case {'posInteger'}
        valueValid = isa(value,'double') && isreal(value) && isscalar(value) ...
            && (value >= 1) && value == floor(value);
        if ~valueValid
            errid = 'globaloptim:typeValueChecker2:NotAPosInteger';
            error(message(errid,propertyName));
        end
    case {'stringsType'}
        strings = valueData;
        valueValid =  ischar(value) && any(strcmp(value,strings));        
        if ~valueValid
            % Format strings for error message
            allstrings = formatCellArrayOfStrings(strings);            
            errid = 'globaloptim:typeValueChecker2:NotAStringsType';
            error(message(errid,propertyName,allstrings));
        end            
    case {'boundedReal'}
        bounds = valueData;
        % Scalar in the bounds
        valueValid = isa(value,'double') && isreal(value) && isscalar(value) ...
            && (value >= bounds(1)) && (value <= bounds(2));
        if ~valueValid
            errid = 'globaloptim:typeValueChecker2:NotABoundedReal';
            error(message(errid,propertyName, sprintf('[%.3g, %.3g]', bounds(1), bounds(2))));
        end        
    case {'functionOrCellArray'}
        % Empty, function handle, string or cell array of functions
        valueValid =  isempty(value) || ischar(value) ...
            || isa(value, 'function_handle') ...
            || (iscell(value) && all(cellfun(@(x) functionChecker(x),value(:))));
        if ~valueValid
            errid = 'globaloptim:typeValueChecker2:notAFunctionOrCellArray';
            error(message(errid,propertyName));
        end
end

%---------------------------------------------------------------------------------
function    allstrings = formatCellArrayOfStrings(strings)
%formatCellArrayOfStrings converts cell array of strings "strings" into an 
% array of strings "allstrings", with correct punctuation and "or" depending
% on how many strings there are, in order to create readable error message.

% To print out the error message beautifully, need to get the commas and "or"s
% in all the correct places while building up the string of possible string values.
allstrings = ['''',strings{1},''''];
for index = 2:(length(strings)-1)
    % add comma and a space after all but the last string
    allstrings = [allstrings, ', ''', strings{index},''''];
end
allstrings = [allstrings,' or ''',strings{end},''''];
    
function isAFun = functionChecker(anInput)
%functionChecker checks whether the input is a string or a function handle.

isAFun = ischar(anInput) || isa(anInput, 'function_handle');