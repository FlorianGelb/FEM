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
            sys = struct('A' , A, 'B' , eye(size(B)), 'C' , C, 'D' , D);
            [rom, info] = ml_ct_ss_bt(sys, opts);
            C = [obj.parameterObj.u0.'];
            c_o = info.W * obj.parameterObj.u0.';
            dt = obj.parameterObj.dt;
            j = 2;
            for i=dt:dt:obj.parameterObj.T 
                c_n = obj.parameterObj.alpha*dt *rom.A* c_o +  dt * rom.B * obj.parameterObj.h(:, j) + c_o;
                c_o = c_n;
                C(:, j) = info.T * c_o;
                j= j +1;
            end
            imagesc(C)
            sol = solution(C, obj.pred.methode, 0, obj.pred);

        end

    end
end