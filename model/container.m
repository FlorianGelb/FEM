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
        
        function sol = euler(obj, A, B, T, c_0)

            C = [T*c_0];
            j = 2;
            dt = obj.parameterObj.dt;
            for i=dt:dt:obj.parameterObj.T 
                c_n = dt *A* c_0 + dt * B*obj.parameterObj.h(:, j) + c_0;
                c_0 = c_n;
                C(:, j) = T*c_0;
                j= j +1;
            end
            sol = solution(C, "", 0, 0); 
        end

    end
end

