% class unit, filename 'unit.m'
%   This class is the abstract parent class for units. 
%   All classes that use 'signal_interface' for input and output transmission need to inherit from the unit class.
% 
% Examples:
%   @code
%   class classA < unit
%   @endcode
% 
% All inheriting subclasses must define the following parameters:
%   nInputs: number of inputs
%   nOutputs: number of outputs
% All inheriting subclasses must implement the function:
%   generate
%
% from @OCG-LAB OPC Framework
%
% 简体中文版注释：
% 类 unit，文件名为 'unit.m'
%   该类是单元的抽象父类。
%   所有使用 'signal_interface' 进行输入和输出传输的类都需要继承自 unit 类。
% 
% 示例：
%   @code
%   class classA < unit
%   @endcode
% 
% 所有继承的子类必须定义以下参数：
%   nInputs：输入数量
%   nOutputs：输出数量
% 所有继承的子类必须实现以下函数：
%   generate
%
% 来自 @OCG-LAB OPC 框架

classdef unit < handle
    properties(Abstract=true, Hidden=true)
        nInputs;    %number of inputs
        nOutputs;   %number of outputs
    end

    properties
        label;      %label of inheriting class
        results;    %for sotring results
    end

    properties(Hidden=true)
        ID;         %ID of inheriting class
    end

    properties(GetAccess=public, SetAccess=protected, Hidden=true)
        inputBuffer = {};
    end

    methods(Abstract)
        varargout = generate(obj, varargin)
    end

    methods
        function obj = unit
            % Constructor function is called when inherited class instantiates an object.
            % Each instantiated object gets a unique ID.
            obj.ID = java.rmi.server.UID();
            % Record the class name.
            obj.label = class(obj);
        end

        function outputs = generateNode(obj)
            if obj.nInputs ~= numel(obj.inputBuffer)
                slog('Number of connected inputs must be equal to the number of module inputs. %s has %d inputs; %d in were given', 'ERR', obj.label, obj.nInputs, numel(obj.inputBuffer));
            end
            outputs{1:obj.nOutputs} = obj.generate(obj.inputBuffer{:});
            obj.inputBuffer = outputs;
        end
    end
end