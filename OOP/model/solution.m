classdef solution

    
    properties
        solution_data;
        methode;
        pod;
        pred;
    end
    
    methods
        function obj = solution(solution_data, methode, pod, pred)
            obj.solution_data = solution_data;
            obj.methode = methode;
            obj.pod = pod;
            obj.pred = pred;
        end
    end
end

