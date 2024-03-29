classdef main
    
    properties
        algorithms = {};
        parameterObj;
        sols = {};
        modes;
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

        function obj = flush_rom(obj)
            indexes = [];
            for i = 1:length(obj.sols)

                if obj.sols{i}.method ~= "Finite Element Method"
                    indexes(end+1) = i;
                end
            end
            obj.sols(indexes) = [];

        
        end
        
        function obj = rom(obj, modes, m)
            %if modes < 10
            %modes = 10;
            %warning("Number of modes has to be at least 10")
            %end
            %if modes > obj.parameterObj.n
            % warning("Number of modes cannot be larger than " + string(obj.parameterObj.n));
            % modes = obj.parameterObj.n;
            %end
            obj.modes = modes;
            
            pod = createPOD(obj.sols{1}, modes, obj.parameterObj);
            obj.sols{end+1} = pod.solve();

            bt = createBALRED(modes, obj.sols{1}, obj.parameterObj);
            obj.sols{end+1} = bt.solve();

            mt = createMODTRUNC(modes, obj.sols{1}, obj.parameterObj);
            obj.sols{end+1} = mt.solve();

            hna = createHNA(modes, obj.sols{1}, obj.parameterObj);
            obj.sols{end+1} = hna.solve();




            
        end


        function obj = solve(obj)
            for i  = 1:length(obj.algorithms)
                solution = obj.algorithms{i}.solve();
                obj.sols{end + 1} = solution;

                if isequal(class(obj.algorithms{i}),'containerFEM')
                    [K M] = obj.algorithms{i}.construct_matrices();
                    im = inv(M);
                    N = obj.parameterObj.alpha * im*K;
                    obj.parameterObj = obj.parameterObj.set_system_matrices(N, eye(size(im)), eye(size(N)),  zeros(size(N, 1), size(im, 2)));

                end
            end
        end
        
        function [E_2] = compare(obj, close_all)
             T = obj.parameterObj.t;
             X = obj.parameterObj.X;
             E_2 = dictionary;
             for i = 1:length(obj.sols)   
                sol = obj.sols{i};
                    
                    if ~isa(sol.pred, "double")
                        %E = sol.solution_data - sol.pred.solution_data;
                        %E_2(sol.method) = norm(E, "fro");
                        sys = ss(obj.parameterObj.A, obj.parameterObj.B, obj.parameterObj.C, obj.parameterObj.D);
                        sys_r = ss(sol.reduced_model.A, sol.reduced_model.B, sol.reduced_model.C, sol.reduced_model.D);
                        sys_e = sys-sys_r;
                        E_2(sol.method) = norm(sys_e, "inf");
%                         E_2 = [];
%                         for t = 1:length(E)
%                             E_2(end + 1) = norm(E(:, t), 2);
%                             %E_max(end+1) = norm(E(:, t), "inf");
%                         end
%                         figure;
%                         plot(T, E_2);
%                         set(gca, 'FontSize', 14);
%                         title("L2 Error of " + sol.method, ' ')
%                         xlabel("Time in s");
%                         ylabel("||Y - Ŷ||_2");
%                        
%                         exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_"+ sol.method + "_"+obj.modes+"_" + obj.parameterObj.n + ".png")
                        %exportgraphics(gcf, "/home/f/Documents/Studienarbeit/images/L2_"+ sol.method + "_"+obj.modes+".png")

                        %figure;
                        %imagesc(T, X, abs(E./sol.solution_data));
                        %set(gca, 'FontSize', 14);
                        %title("Relative Error" + sol.method, ' ')
                        %xlabel("Time in s");
                        %ylabel("Length in cm");
                        %colormap turbo;
                        %colorbar;
                        %caxis([0, 1]); 
                        
                        %exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/abs_"+ sol.method + "_"+obj.modes+ "_" + obj.parameterObj.n + ".png")
                        %exportgraphics(gcf, "/home/f/Documents/Studienarbeit/images/abs_"+ sol.method + "_"+obj.modes+".png")
                        
                        %figure;
                        %plot(T,  E_max);
                        %set(gca, 'FontSize', 14);
                        %title("Max Error for " + sol.method, " ")
                        %xlabel("Time in s");
                        %ylabel("||Y - Ŷ||_{max}");
                        
                        %exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/max_"+ sol.method + "_"+obj.modes+ "_" + obj.parameterObj.n +".png")
                        %exportgraphics(gcf, "/home/f/Documents/Studienarbeit/images/max_"+ sol.method + "_"+obj.modes+".png")
                        
                        if close_all
                            close all
                        end
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
                        set(gca, 'FontSize', 18);
                    end
                colorbar;
                caxis([min_, max_]);  
            end
        end
    end
end

