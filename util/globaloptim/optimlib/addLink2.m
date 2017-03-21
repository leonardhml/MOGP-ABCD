function [taggedString,openTag,closeTag] = addLink2(linkedText, ... 
                        toolboxName,mapFileName,linkAnchor,isCSH)
%

%addLink2 add a hyperlink to a string for display in the MATLAB Command
%Window.
%
%   taggedString = addLink2(linkedText,linkDestination) takes an input
%   string (linkedText) and wraps it in html tags that execute a MATLAB
%   command to open the documentation browser to a specified location
%   (linkDestination) in the Optimization Toolbox documentation. The result
%   (taggedString) can be inserted in any text printed to the MATLAB
%   Command Window (e.g. error, MException, warning, fprintf).
%
%   taggedString = addLink2(linkedText,linkDestination,toolboxName) directs
%   the MATLAB command to open the documentation for the specified toolbox.

%   Copyright 2009-2015 The MathWorks, Inc.

if feature('hotlinks') && ~isdeployed;      
    % Assemble helpview arguments
    helpviewPath = sprintf('''%s/%s/%s''',docroot,toolboxName,mapFileName);

    windowType = '';
    if isCSH
        windowType = ',''CSHelpWindow''';
    end
    % Create explicit char array so as to avoid translation
    openTag = sprintf('<a href = "matlab: helpview(%s,''%s''%s);">',...
        helpviewPath,linkAnchor,windowType);
    closeTag = '</a>';
    taggedString = [openTag linkedText closeTag];
else
    taggedString = linkedText;
    openTag = '';
    closeTag = '';
end
