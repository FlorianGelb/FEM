classdef createPOD < container
    
    properties
        snapshots;
        energie;
    end
    
    methods
        function obj = createPOD(snapshots,modes, parameterObj)
            obj = obj@container(parameterObj);
            obj.snapshots = snapshots;
            obj.energie = modes;   
        end
        
        function sol = solve(obj)
            snapshots = obj.snapshots.solution_data;
            [U S V] = svd(snapshots, "econ");
            E = diag(S)/trace(S);
            %figure;
            %plot(100*E, "*", 'MarkerSize',15);
            %title("Variance Captured by Mode");
            %xlabel("Mode Number");
            %ylabel("Variance in %");
            %set(gca, 'FontSize', 25);

            %n_mode = 0;
            %cum_E = 0;
            %for i = 1:length(E)
            %    cum_E = cum_E + E(i);
            %    n_mode = n_mode + 1;
            %    if cum_E >= obj.energie
            %    break
            %    end
            %end
            n_mode = obj.energie;
            phi = U(:, 1:n_mode);
            phi_xx = [];
            for i = 1:n_mode
                phi_xx(:,i) = obj.spectral_der(phi(:, i));
            end

            a0 = phi.' * snapshots(:, 1);
            dt = obj.parameterObj.dt;
            [t, asol] = ode45("a_rhs", obj.parameterObj.t, a0,[], phi, phi_xx, obj.parameterObj.alpha, dt, obj.parameterObj.h);
            
            usol = asol*phi.'
            rom = struct("A", obj.parameterObj.alpha*phi.'*phi_xx, "B", phi.', "C", phi, "D", obj.parameterObj.D);
            sol = solution(usol.', "Proper Orthorgonal Decomposition", obj.snapshots, rom);
        end

        function der = spectral_der(obj, v)
        der = ifft(-1 * obj.parameterObj.k.^2.* fft(v));
        end
    end
end

