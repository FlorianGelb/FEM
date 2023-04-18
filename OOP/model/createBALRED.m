classdef createBALRED < container
    properties
    modes;
    pred;
    end
    methods
        function obj = createBALRED(modes, pred, parameterObj)
             obj = obj@container(parameterObj);
             obj.modes = modes;
             obj.pred = pred;
        end
        
        function sol = solve(obj)
            A = obj.parameterObj.A;
            B = obj.parameterObj.B;
            C = obj.parameterObj.C;
            D = obj.parameterObj.D;
            opts = ml_morlabopts('ml_ct_ss_bt');
            opts.StoreProjection = 1;
            opts.OrderComputation = 'Order';
            opts.Order            = obj.modes;
            sys = struct('A' , A, 'B' , B, 'C' , C, 'D' , D);
            [rom, info] = ml_ct_ss_bt(sys, opts);
            c_0 = info.W * obj.parameterObj.u0.';
            sol = obj.euler(rom.A, rom.B, info.T, c_0);
            sol.method = "Balanced Reduction";
            sol.pred = obj.pred;
            sol.reduced_model = rom;
        end

    end
end