classdef solution

    
    properties
        solution_data;
        methode;
        pod;
    end
    
    methods
        function obj = solution(solution_data, methode, pod)
            obj.solution_data = solution_data;
            obj.methode = methode;
            obj.pod = pod;
        end
    end
end

