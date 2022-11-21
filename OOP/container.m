classdef container
    
    properties
        parameterObj;
    end
    
    methods
        function obj = container(parameterObj)
            obj.parameterObj = parameterObj;
        end
        
        function sol = solve(obj)
            sol = NaN; 
        end

    end
end

