classdef createPOD
    
    properties
        snapshots;
        energie;
        parameterObj;
    end
    
    methods
        function obj = createPOD(snapshots,energie, parameterObj)
            obj.snapshots = snapshots;
            obj.energie = energie;
            obj.parameterObj = parameterObj;      
        end
        
        function sol = solve(obj)
            snapshots = obj.snapshots.solution_data;
            [U S V] = svd(snapshots, "econ");
            E = diag(S)/trace(S);
            plot(E)
            n_mode = 0;
            cum_E = 0;
            for i = 1:length(E)
                cum_E = cum_E + E(i);
                n_mode = n_mode + 1;
                if cum_E >= obj.energie
                break
                end
            end
            
            phi = U(:, 1:n_mode);
            phi_xx = [];
            for i = 1:n_mode
                phi_xx(:,i) = obj.spectral_der(phi(:, i));
            end

            a0 = phi.' * snapshots(:, 1);
            [t, asol] = ode45("a_rhs", obj.parameterObj.t, a0,[], phi, phi_xx, obj.parameterObj.alpha);
            
            usol = zeros(size(phi, 1), length(obj.parameterObj.t)).';           
            for j = 1:length(t)
                for i = 1:n_mode
                    usol(j,:) = usol(j, :) + asol(j,i)*phi(:,i).';
                end
            end
            sol = solution(usol.', obj.snapshots.methode, 1, obj.snapshots);
        end

        function der = spectral_der(obj, v)
        der = ifft(-1 * obj.parameterObj.k.^2.* fft(v));
        end
    end
end

