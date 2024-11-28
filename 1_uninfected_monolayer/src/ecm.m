classdef ecm

    properties %(Access = public)
        ecm_size = [0 0 0] % ECM Size µm
    end

    methods
        function obj = setECMSize(obj,x,y,z)
            % Value class method that modifies the object must return the
            % modified object
            arguments
                obj 
                x   {mustBeNumeric}
                y   {mustBeNumeric}
                z   {mustBeNumeric}
            end

            if ~isa(obj, 'ecm')
                error('Error: the input is not an ecm object');
            end

            obj.ecm_size = [x y z];
        end
        
    end




end