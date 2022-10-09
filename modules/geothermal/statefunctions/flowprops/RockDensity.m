classdef RockDensity < StateFunction

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RockDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures', 'Temperature'});
            gp.label = '\rho_R';
        end
        
        %-----------------------------------------------------------------%
        function rhoR = evaluateOnDomain(prop, model, state)
            [p, T] = prop.getEvaluatedDependencies(state, 'PhasePressures', ...
                                                          'Temperature'   );
            rhoR = prop.evaluateFunctionOnDomainWithArguments(model.rock.rhoR, p, T);
        end
    end
    
end