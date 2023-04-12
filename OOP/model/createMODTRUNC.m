classdef createMODTRUNC < container
    properties
    modes;
    pred;
    end
    methods
        function obj = createMODTRUNC(modes, pred, parameterObj)
             obj = obj@container(parameterObj);
             obj.modes = modes;
             obj.pred = pred;
        end

        
        function sol = solve(obj)
            A = obj.parameterObj.A;
            B = obj.parameterObj.B;
            C = obj.parameterObj.C;
            D = obj.parameterObj.D;

            eigen_values = sort(real(eig(A)))
            size(eigen_values)
            alpha = (eigen_values(end-(obj.modes)) + eigen_values(end-(obj.modes-1)))/2
            size(alpha)
            
            opts = ml_morlabopts('ml_ct_ss_mt');
            opts.StoreProjection = 1;
            opts.OrderComputation = 'Alpha';
            opts.Alpha = alpha;
            sys = struct('A' , A, 'B' , B, 'C' , C, 'D' , D);
            [rom, info] = ml_ct_ss_mt(sys, opts)
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
            sol = solution(C, obj.pred.methode, 0, obj.pred);

        end

    end
end