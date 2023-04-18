classdef main
    
    properties
        algorithms = {};
        parameterObj;
        sols = {};
    end
    
    methods
        function obj = main(algorithms, L, u0, h, T, n, nt, alpha)
            obj.parameterObj = parameters(L, u0, h, T, n, nt, alpha);
            for i  = 1:length(algorithms)
                if algorithms(i) == "FT"
                    obj.algorithms{end + 1} = containerFT(obj.parameterObj);
                elseif algorithms(i) == "FEM"
                    obj.algorithms{end + 1} = containerFEM(obj.parameterObj);
                end
            end
        end
        
        function obj = rom(obj, energie, modes)
            if modes < 10
            modes = 10;
            warning("Number of modes has to be at least 10")
            end
            if modes > obj.parameterObj.n
             warning("Number of modes cannot be larger than " + string(obj.parameterObj.n));
             modes = obj.parameterObj.n;
            end
            for i = 1:length(obj.sols)
              % pod = createPOD(obj.sols{i}, energie, obj.parameterObj);
              
               bt = createBALRED(modes, obj.sols{i}, obj.parameterObj);
               mt = createMODTRUNC(modes, obj.sols{i}, obj.parameterObj);
               hna = createHNA(modes, obj.sols{i}, obj.parameterObj);
               %obj.sols{end+1} = pod.solve();
               obj.sols{end+1} = bt.solve();
               obj.sols{end+1} = mt.solve();
               obj.sols{end+1} = hna.solve();
            end
        end


        function obj = solve(obj)
            for i  = 1:length(obj.algorithms)
                solution = obj.algorithms{i}.solve();
                obj.sols{end + 1} = solution;

                if isequal(class(obj.algorithms{i}),'containerFEM')
                    [K M] = obj.algorithms{i}.construct_matrices();
                    im = inv(M);
                    N = im*K;
                    obj.parameterObj = obj.parameterObj.set_system_matrices(N, eye(size(im)), eye(size(N)),  zeros(size(N, 1), size(im, 2)));

                end
            end
        end
        
        function compare(obj)
             for i = 1:length(obj.sols)   
                sol = obj.sols{i};
                    figure(13)
                    if ~isa(sol.pred, "double")
                        E = sol.solution_data - sol.pred.solution_data;
                        E_2 = []
                        E_max = [];
                        for t = 1:length(E)
                            E_2(end + 1) = norm(E(:, t), 2);
                            E_max(end+1) = norm(E(:, t), "inf");
                        end
                        plot(E_2);
                        figure(14);
                        imagesc(E);
                        figure(15);
                        plot(E_max);
                    end 
             end
        end

        function plot(obj)  
            T = obj.parameterObj.t;
            X = obj.parameterObj.X;
            min_ = Inf;
            max_ = -Inf;
            
            
            
            for i = 1:length(obj.sols)
                     
                    sol = obj.sols{i};
                    if isstruct(sol.reduced_model)
                        figure(100+i);
                        sys = struct('A' , obj.parameterObj.A, 'B' ,  obj.parameterObj.B, 'C' ,  obj.parameterObj.C, 'D' ,  obj.parameterObj.D);
                        ml_frobeniusplot(sol.reduced_model, sys);
                        ml_sigmaplot(sol.reduced_model, sys, struct('DiffMode', 'abs'));
                    end
                    buffer = {};
                    if isa(sol.pred, "double")
                        continue;
                    end
                    while ~isa(sol.pred, "double")
                        buffer{end + 1} = obj.sols{i};
                        sol = sol.pred;
                        if isa(sol.pred, "double")
                            buffer{end + 1} = sol;
                            break;
                        end
                        
                    end
                    
           
                    figure(i)
                    for j = 1:length(buffer)
                        subplot(1, length(buffer) , j);
                        h = imagesc(T, X,buffer{j}.solution_data);
                        mx = max(max(buffer{j}.solution_data));
                        if mx > max_
                            max_ = mx;
                        end
                        mn = min(min(buffer{j}.solution_data));
                        if mn < min_
                            min_ = mn;
                        end
                
                        title(buffer{j}.method);
                        colormap turbo;        
                    end
                colorbar;
                caxis([min_, max_]);  
            end
        end
    end
end

