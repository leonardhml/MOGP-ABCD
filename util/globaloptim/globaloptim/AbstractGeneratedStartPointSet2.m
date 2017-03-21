classdef (Hidden) AbstractGeneratedStartPointSet2 < AbstractStartPointSet2
%AbstractGeneratedStartPointSet2 Generated start point set class.
%   AbstractGeneratedStartPointSet2 is an abstract class representing a
%   generated start point set. You cannot create instances of this class 
%   directly.  You must create a derived class such as RandomStartPointSet2,
%   by calling the class constructor.
%
%   All start point sets that inherit from AbstractGeneratedStartPointSet2
%   have the following properties and method:
%
%   AbstractGeneratedStartPointSet2 properties:
%       NumStartPoints      - number of start points in the set
%       ArtificialBound     - artificial bound for unbounded problems
%
%   AbstractGeneratedStartPointSet2 method:
%       list                - lists the start points in the set
%
%   See also RANDOMSTARTPOINTSET2, ABSTRACTSTARTPOINTSET2, CUSTOMSTARTPOINTSET2

%   Copyright 2009-2011 The MathWorks, Inc.

    properties
%NUMSTARTPOINTS Number of start points in the set
%   The NumStartPoints property indicates the number of the start
%   points. The default number of start points is 10.
%
%   See also ABSTRACTGENERATEDSTARTPOINTSET2, LIST
        NumStartPoints = 10;
        
%ARTIFICIALBOUND Artificial bound for unbounded problems
%   The ArtificialBound property sets an arbitrary bound for the variables
%   in the problem only if there is no specified upper or lower bound. A
%   GeneratedStartPointSet will create start points in a bounded box. If
%   the upper (lower) bound is not specified for a variable, the
%   GeneratedStartPointSet sets the upper (lower) bound to be
%   ArtificialBound (-ArtificialBound) for that variable. 
%
%   See also ABSTRACTGENERATEDSTARTPOINTSET2, LIST
        ArtificialBound = 1000;        
    end
    methods
        % Constructor
        function obj = AbstractGeneratedStartPointSet2(varargin)
            % First argument can be a RandomStartPointSet2, or parameter of
            % a parameter-value pair.
            firstInputObj = false;
            if nargin > 0
                value = varargin{1};
                if isa(value,'AbstractGeneratedStartPointSet2')
                    if strcmp(class(value),class(obj))
                        obj = value;
                    else
                        obj.NumStartPoints = value.NumStartPoints;
                        obj.ArtificialBound = value.ArtificialBound;
                    end
                    firstInputObj = true;
                end
                if nargin > 1 && firstInputObj
                    agtpsProp = varargin(2:end);
                elseif nargin == 1 && ~firstInputObj
                    errid = 'globaloptim:AbstractGeneratedStartPointSet2:invalidInput';
                    error(message(errid));
                else
                    agtpsProp = varargin;
                end
                if (nargin > 1) && (mod(length(agtpsProp),2) ~= 0)
                    error(message('globaloptim:AbstractGeneratedStartPointSet2:notAPair'));
                end
                while length(agtpsProp) >=2
                    prop = agtpsProp{1};
                    val = agtpsProp{2};
                    agtpsProp = agtpsProp(3:end);
                    if ~ischar(prop)
                       error(message('globaloptim:AbstractGeneratedStartPointSet2:incorrectPropertyEntry'));
                    end
                    PName = i_isPublicExact(obj, prop, 'set');
                    if ~isempty(PName)
                        obj.(PName) = val;
                    else
                        % No match at all
                        error(message('globaloptim:AbstractGeneratedStartPointSet2:invalidPropertyName', prop));
                    end
                end
            end
        end
        function obj = set.NumStartPoints(obj,value)
            % error check for NumStartPoints property
            typeValueChecker2('posInteger',value,'NumStartPoints');
            obj.NumStartPoints = value;
        end
        function obj = set.ArtificialBound(obj,value)
            % error check for ArtificialBound property
            typeValueChecker2('posReal',value,'ArtificialBound');
            obj.ArtificialBound = value;
        end
    end
    methods (Abstract, Access = protected)
        startPoints = generate(obj,problem)
    end      
end

function CorrectIdent = i_isPublicExact(obj, Ident, Type)
% Check whether an identifier matches exactly to a public property or
% method.

publicIdents = i_getPublic(obj,Type);
CorrectIdent = '';
if strmatch(Ident,publicIdents,'exact');
    CorrectIdent = Ident;
end
end

function Names = i_getPublic(obj, Type)
% Return the list of public names for the specified type of operation.
% Type may be one of 'get', 'set' or 'call'.

% The list of public identifiers will depend on the exact class of the
% object, so a list is held for each unique class.

persistent PublicIdents
if isempty(PublicIdents)
    % Initialise the database of lists-per-class to be empty
    PublicIdents = struct;
end

Cname = class(obj);
if ~isfield(PublicIdents, Cname)
    % Initialise the lists of public identifiers for this class
    C = metaclass(obj);
    P = [C.Properties{:}];
    PropGet = strcmp({P.GetAccess}, 'public');
    PropSet = strcmp({P.SetAccess}, 'public');
    PropNames = {P.Name};

    M = [C.Methods{:}];
    MethAccess = strcmp({M.Access}, 'public');
    MethNames = {M.Name};

    PublicIdents.(Cname).get = PropNames(PropGet);
    PublicIdents.(Cname).set = PropNames(PropSet);
    PublicIdents.(Cname).call = MethNames(MethAccess);
end

Names = PublicIdents.(Cname).(Type);
end
