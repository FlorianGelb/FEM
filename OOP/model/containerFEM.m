classdef containerFEM < container

    
    properties
        delta_nodes;
        nodes;       
    end
    methods
        function obj = containerFEM(parameterObj)
             obj = obj@container(parameterObj);
             obj.delta_nodes = obj.parameterObj.dx;
             obj.nodes = obj.parameterObj.n;
        end

        function sol = solve(obj)
            nodes = obj.nodes;
            [K M] = obj.construct_matrices();
            N = inv(M) * K;
            sol = obj.euler(N, eye(size(M)), eye(size(M)), obj.parameterObj.u0.');
            sol.method = "Finite Element Method";
        end
    
        function [K M] = construct_matrices(obj)
             delta_nodes = obj.delta_nodes;
             nodes = obj.nodes;
             ii = (2/3) * delta_nodes;
             ij = (1/6) * delta_nodes;
             
             K = zeros(nodes);
             M = zeros(nodes);

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

