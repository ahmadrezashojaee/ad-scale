classdef PhaseThermalEnergy < StateFunction
       
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = PhaseThermalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseMass', 'PhaseInternalEnergy'});
            gp.label = 'E_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function energy = evaluateOnDomain(prop,model, state)
            [mass, u] = prop.getEvaluatedDependencies(state, 'PhaseMass'          , ...
                                                             'PhaseInternalEnergy'); 
            nph    = model.getNumberOfPhases();
            energy = cell(1, nph);
            for i = 1:nph
                energy{i} = mass{i}.*u{i};
            end
        end       
    end
end

