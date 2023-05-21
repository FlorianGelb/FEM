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

        
        function sol = solve(obj, time)
            A = obj.parameterObj.A;
            B = obj.parameterObj.B;
            C = obj.parameterObj.C;
            D = obj.parameterObj.D;

            eigen_values = sort(real(eig(A)));


            if obj.modes == obj.parameterObj.n
            alpha = eigen_values(1) - 10;
            else
                alpha = (eigen_values(end-(obj.modes)) + eigen_values(end-(obj.modes-1)))/2;
            end

            
            opts = ml_morlabopts('ml_ct_ss_mt');
            opts.StoreProjection = 1;
            opts.OrderComputation = 'Alpha';
            opts.Alpha = alpha;
            sys = struct('A' , A, 'B' , B, 'C' , C, 'D' , D);
            rom_time = [];
            
            if time
                for i = 1:10
                    f = @() ml_ct_ss_mt(sys, opts);
                    rom_time(end+1) = timeit(f);
                end

                sol = solution(NaN, "Modal Truncation", obj.pred, sys);
                sol.rom_time = rom_time;
            else
                [rom, info] = ml_ct_ss_mt(sys, opts);
                c_0 = info.W * obj.parameterObj.u0.';
                sol = obj.euler(rom.A, rom.B, info.T, c_0);
                sol.method = "Modal Truncation";
                sol.pred = obj.pred;
                sol.reduced_model = rom;
                sol.rom_time = rom_time;
            end

        end

    end
end