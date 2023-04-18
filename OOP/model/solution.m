classdef solution

    
    properties
        solution_data;
        method;
        pred;
        reduced_model;
    end
    
    methods
        function obj = solution(solution_data, method, pred, mred)
            obj.solution_data = solution_data;
            obj.method = method;
            obj.pred = pred;
            obj.reduced_model = mred;
        end
    end
end

