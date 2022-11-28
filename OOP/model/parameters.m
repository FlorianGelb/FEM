classdef parameters   
    properties
        L;
        u0;
        T;
        n;
        dx;
        dt;
        X;
        t;
        alpha;
        k;
        nt;
        n_nodes;
    end
    
    methods
        function obj = parameters(L, u0, T, n, nt, alpha)
            if L <= 0 || T <= 0 || n <= 0 || nt <= 0 || alpha < 0 
                E = MException("parameters:WrongParameter", "Invalid parameter set");
                throw(E);
            end

            obj.L = L;
            obj.u0 = u0;
            obj.T = T;
            obj.n = n;
            obj.nt = nt;
            obj.dx = L/(n-1);
            obj.alpha = alpha;
            obj.X = linspace(0, L, n);
            obj.t = linspace(0, T, nt);
            obj.k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
        end    

    end
end

