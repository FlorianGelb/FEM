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
        
        function sol = solve(obj, time)
            snapshots = obj.snapshots.solution_data;
            n_mode = obj.energie;
            rom_time = [];
            if time
                for i = 1:10
                    tic;
                    [U S V] = svd(snapshots, "econ");
                    phi = U(:, 1:n_mode);
                    phi_xx = [];
                    for i = 1:n_mode
                        phi_xx(:,i) = obj.spectral_der(phi(:, i));
                    end
                    rom_time(end+1) = toc;
                end
                sol = solution(NaN, "Proper Orthogonal Decomposition", obj.snapshots, NaN);
                sol.rom_time = rom_time;
            else
            
            [U S V] = svd(snapshots, "econ");
            phi = U(:, 1:n_mode);
            phi_xx = [];
            for i = 1:n_mode
                phi_xx(:,i) = obj.spectral_der(phi(:, i));
            end

            a0 = phi.' * snapshots(:, 1);
            dt = obj.parameterObj.dt;
            sol = obj.euler(obj.parameterObj.alpha*phi.'*phi_xx, phi.', phi, a0);
                
            rom = struct("A", obj.parameterObj.alpha*phi.'*phi_xx, "B", phi.', "C", phi, "D", obj.parameterObj.D);
            sol.method = "Proper Orthogonal Decomposition";
            sol.pred = obj.snapshots;
            sol.reduced_model = rom;
            sol.rom_time = rom_time;
            end
        end

        function der = spectral_der(obj, v)
        der = ifft(-1 * obj.parameterObj.k.^2.* fft(v));
        end
    end
end

