classdef createHNA < container
    properties
    modes;
    pred;
    end
    methods
        function obj = createHNA(modes, pred, parameterObj)
             obj = obj@container(parameterObj);
             obj.modes = modes;
             obj.pred = pred;
        end

        
        function sol = solve(obj)
            A = obj.parameterObj.A;
            B = obj.parameterObj.B;
            C = obj.parameterObj.C;
            D = obj.parameterObj.D;
            
            opts = ml_morlabopts('ml_ct_ss_hna');
            opts.StoreProjection = 1;
            opts.OrderComputation = 'Order';
            opts.Order            = obj.modes;
            sys = struct('A' , A, 'B' , B, 'C' , C, 'D' , D);
            [rom, info] = ml_ct_ss_hna(sys, opts)
            c_o = zeros(obj.modes, 1);
            C = [c_o];
            dt = obj.parameterObj.dt;
            j = 2;
            
            size(ss2tf(rom.A, rom.B, rom.C, rom.D, ))

            for i=dt:dt:obj.parameterObj.T 
                c_n = obj.parameterObj.alpha*dt *rom.A* c_o +  dt * rom.B * obj.parameterObj.h(:, j) + c_o;
                c_o = c_n;
                C(:, j) = c_o;
                j= j +1;
            end
            sol = solution(C, obj.pred.methode, 0, obj.pred);

        end

    end
end