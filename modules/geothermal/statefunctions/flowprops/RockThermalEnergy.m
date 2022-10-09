classdef RockThermalEnergy < StateFunction

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RockThermalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'RockMass', 'RockInternalEnergy'});
            gp.label = 'E_R';
        end
        
        %-----------------------------------------------------------------%
        function uR = evaluateOnDomain(prop, model, state)
            [massR, uR] = prop.getEvaluatedDependencies(state, 'RockMass'          , ...
                                                               'RockInternalEnergy');
            uR = massR.*uR;
        end 
    end
    
end