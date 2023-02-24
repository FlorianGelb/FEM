classdef containerFEM < container

    
    properties
        delta_nodes;
        nodes;
       
    end

    methods(Static)
        function val = triag(x, delta_x, x0, x1, x2)
            if x0 <= x && x < x1
                val = (x-x0)/delta_x;
            elseif x1 <= x && x < x2
                val = (x2 - x) / delta_x;
            else
                val = 0;
            end
        end
    end
    methods
        function obj = containerFEM(parameterObj)
             obj = obj@container(parameterObj);
             obj.delta_nodes = 2 * obj.parameterObj.dx;
             obj.nodes = floor(obj.parameterObj.L / obj.delta_nodes + 1);
        end

        function sol = solve(obj)
            nodes = obj.nodes;
            [F K M] = obj.construct_matrices();
            im = inv(M);
            %im(1,:) = [0, zeros(1, nodes-1)];
            %im(nodes,:) = [zeros(1, nodes-1), 0];
            t = obj.parameterObj.t;
            C = zeros(nodes, uint16(obj.parameterObj.nt));
            dt = obj.parameterObj.dt;
            c_o = im*F;
            C(:, 1) = c_o;
            j = 2;
            N = im * K;

            for i=dt:dt:obj.parameterObj.T 
                c_n = obj.parameterObj.alpha*dt *N* c_o + dt * im * obj.generate_h_disc(i) + c_o;
                c_o = c_n;
                C(:, j) = c_o;
                j= j +1;
            end
            S = [];
            
            for t = 1:1:obj.parameterObj.nt
                c = C(:, t);
                interp_c = interp1(linspace(0, obj.parameterObj.L, obj.nodes), c, obj.parameterObj.X);
                S(:,t) = interp_c;
            end
            S
            sol = solution(S, "FEM", 0, 0);
        end
        
        function t = generate_triags(obj, i)
            delta_nodes = obj.delta_nodes;
            t = [];
            for j = 1:length(obj.parameterObj.X) 
                t(end+1) = containerFEM.triag(obj.parameterObj.X(j), delta_nodes, (i-1)*delta_nodes, (i)*delta_nodes, (i+1)*delta_nodes);
            end
        end


        function V = generate_h_disc(obj, j)
             dt = obj.parameterObj.dt;
             l = int32(j / dt);
             nodes = obj.nodes;
             V = zeros(nodes, 1);
             h = obj.parameterObj.h(:, l).';
             delta_nodes = obj.delta_nodes;
             for i = 1:nodes
                  t = obj.generate_triags(i);
                  V(i) =  trapz(obj.parameterObj.X,t.*h);
             end


        end
        
        function [F K M] = construct_matrices(obj)
             delta_nodes = obj.delta_nodes;
             nodes = obj.nodes;
             ii = (2/3) * delta_nodes;
             ij = (1/6) * delta_nodes;
             F = zeros(nodes, 1);
             f = obj.parameterObj.u0;
             
             K = zeros(nodes);
             M = zeros(nodes);

            for i = 1:nodes
              t = obj.generate_triags(i);
              size(t)
              size(f)
              F(i) =  trapz(obj.parameterObj.X,t.*f);
            end

            figure(499)
            plot(f)
            for i = 1:nodes
                K(i,i) = -2/delta_nodes;
                M(i, i) = ii;
                if i-1 > 0
                    K(i, i-1) = 1/delta_nodes;
                    M(i, i-1) =  ij;
                end
                if i + 1 < nodes+1
                    K(i, i+1) = 1/delta_nodes;
                    M(i, i+1) =   ij;
                end
            end


            

        end


    end
end

