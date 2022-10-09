classdef PhaseMass < StateFunction
       
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = PhaseMass(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'Density', 'PoreVolume'});
            gp.label = 'M_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function mass = evaluateOnDomain(prop,model, state)
            [rho, pv] = prop.getEvaluatedDependencies(state, 'Density'   , ...
                                                             'PoreVolume');
            nph  = model.getNumberOfPhases();
            mass = cell(1, nph);
            for i = 1:nph
                mass{i} = rho{i}.*pv;
            end
        end       
    end
end
