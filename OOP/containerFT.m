classdef containerFT < container    
    properties
        
    end
    
    methods
        function obj = containerFT(parameterObj)
            obj = obj@container(parameterObj);
        end
        
        function sol = solve(obj)
            t = obj.parameterObj.t;
            [t, fft_sol] = ode45(@(t,fft_sol)rhs_1dHeatFT(t,obj.parameterObj.k, obj.parameterObj.alpha, fft_sol),t,fft(obj.parameterObj.u0));
 

            usol = [];
            for j = 1:length(t)
                usol(j, :) = ifft(fft_sol(j,:));    
            end
            sol = solution(usol.', "FT", 0);
        end
    end
end

