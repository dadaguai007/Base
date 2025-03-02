classdef module < unit

    properties(Abstract)
        nInputs;
        nOutputs;
    end
    
    properties
        interunits = {};
        generatingOrder = [];
    end

    properties(SetAccess=protected, GetAccess=protected)
        outputBuffer;
    end

    methods
        function obj = module
        end

        function varargout = generate(obj, varargin)
            obj.interunits{obj.generatingOrder(1)}.inputBuffer = varargin;
            obj.outputBuffer = generateNode(obj.interunits{obj.generatingOrder(1)});
            for idx = obj.generatingOrder(2:end)
                obj.interunits{idx}.inputBuffer = obj.outputBuffer;
                obj.outputBuffer = generateNode(obj.interunits{idx});
            end
            varargout{1:obj.nOutputs} = obj.outputBuffer;
            if obj.nOutputs == 1
                varargout = varargout{:};
            end
        end
    end

end