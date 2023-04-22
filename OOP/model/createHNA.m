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
            [rom, info] = ml_ct_ss_hna(sys, opts);  
            sol = obj.euler(rom.A, rom.B, rom.C, rom.B*obj.parameterObj.u0.');
            sol.method =  "Hankel Norm Approximation";
            sol.reduced_model = rom;
            sol.pred = obj.pred;
        end

    end
end