function classes = findSubClasses2(packageName, superclassName)
%findSubClasses2   Find sub-classes within a package
%
%   CLASSES = findSubClasses2(PACKAGE, SUPERCLASS) is an array of meta.class
%   objects, each element being a sub-class of SUPERCLASS and a member of
%   the given PACKAGE.
%
%   Note that only non-abstract classes are returned.
%
%   Example
%      classes = optim.internal.findSubClasses2('optim.options', 'optim.options.SolverOptions2')

%   Copyright 2015 The MathWorks, Inc.

narginchk( 2, 2 ) ;

% Get the package object
package = meta.package.fromName(packageName);
if isempty(package)
    warning('optimlib:internal:findSubClasses2', 'Package not found.');
    classes = {};
else
    % For each class in the package ...
    %  1. check for given super-class
    %  2. check for abstract classes
    classes = package.ClassList;
    keep = arrayfun( ...
        @(x) isAClass(superclassName, x.SuperClasses ) && ~x.Abstract, classes);
    
    % Return list of non-abstract classes that sub-class the given super-class
    classes = classes(keep);
end

function tf = isAClass( className, list )
% Check the LIST of classes and their superclasses for given CLASSNAME
tf = false;
for i = 1:length( list )
    tf = strcmp( className, list{i}.Name ) || isAClass( className, list{i}.SuperClasses );
    if tf
        break
    end
end