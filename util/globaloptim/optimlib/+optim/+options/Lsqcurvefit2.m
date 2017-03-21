classdef (Sealed) lsqcurvefit2 < optim.options.lsqncommon2
%

%lsqcurvefit2 Options for lsqcurvefit2
%
%   The OPTIM.OPTIONS.lsqcurvefit2 class allows the user to create a set of
%   options for the lsqcurvefit2 solver. For a list of options that can be set,
%   see the documentation for lsqcurvefit2.
%
%   OPTS = OPTIM.OPTIONS.lsqcurvefit2 creates a set of options for lsqcurvefit2
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.lsqcurvefit2(PARAM, VAL, ...) creates a set of options
%   for lsqcurvefit2 with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.lsqcurvefit2(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MultiAlgorithm2, OPTIM.OPTIONS.SolverOptions2
    
%   Copyright 2012-2015 The MathWorks, Inc.    
       
    properties (Hidden)
        SolverName = 'lsqcurvefit2';
    end       

    properties (Hidden, SetAccess = private, GetAccess = public)
        
        % New version property added in third version
        LsqcurvefitVersion
    end
    
    methods (Hidden)
        
        function obj = lsqcurvefit2(varargin)
%lsqcurvefit2 Options for lsqcurvefit2
%
%   The OPTIM.OPTIONS.lsqcurvefit2 class allows the user to create a set of
%   options for the lsqcurvefit2 solver. For a list of options that can be set,
%   see the documentation for lsqcurvefit2.
%
%   OPTS = OPTIM.OPTIONS.lsqcurvefit2 creates a set of options for lsqcurvefit2
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.lsqcurvefit2(PARAM, VAL, ...) creates a set of options
%   for lsqcurvefit2 with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.lsqcurvefit2(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MultiAlgorithm2, OPTIM.OPTIONS.SolverOptions2
            
            % Call the superclass constructor
            obj = obj@optim.options.lsqncommon2(varargin{:});
               
            % Record the class version; Update property 'LsqcurvefitVersion'
            % instead of superclass property 'Version'.
            obj.Version = 2;
            obj.LsqcurvefitVersion = 4;    
        end
        
    end
    
    methods (Static)
        
        function obj = loadobj(obj)
    
            % Objects saved in R2013a will come in as structures. 
            if isstruct(obj) && obj.Version == 1

                % Save the existing structure
                s = obj;
                
                % Create a new object
                obj = optim.options.lsqcurvefit2;
                
                % Call the superclass method to upgrade the object
                obj = upgradeFrom13a(obj, s); 
                
                % The SolverVersion property was not present in 13a. We
                % clear it here and the remainer of loadobj will set it
                % correctly.
                obj.LsqcurvefitVersion = [];
                
            end
            
            % Upgrading to 13b
            % Update Levenberg-Marquardt and finite difference options
            if obj.Version < 2
                obj = upgradeLMandFinDiffTo13b(obj);
            end

            % Set the superclass version number
            obj.Version = 2;

            % Upgrading to 15b            
            % Introduce LsqcurvefitVersion field
            if isempty(obj.LsqcurvefitVersion)
                % Use 'LsqcurvefitVersion' property instead of 'Version'
                % property because 'Version' is for the superclass and
                % 'LsqcurvefitVersion' is for this (derived) class. However,
                % 'LsqcurvefitVersion' was added in the second version, we
                % check only for the second version and add this property.
                % For all other version, check only the 'LsqcurvefitVersion'
                % property.
                obj.LsqcurvefitVersion = 2; % update object
            end
            
            % Upgrading to 16a
            % Add TolFunValue field
            if obj.LsqcurvefitVersion < 4
                obj = upgradeTo16a(obj, 'optim.options.lsqcurvefit2');
            end
            
            % Set the version number
            obj.LsqcurvefitVersion = 4;            
            
        end
        
    end

end
    
