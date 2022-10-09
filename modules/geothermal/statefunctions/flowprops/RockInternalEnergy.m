classdef RockInternalEnergy < StateFunction

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RockInternalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('Temperature');
            gp = gp.dependsOn({'pressure'}, 'state');
            gp.label = 'u_R';
        end
        
        %-----------------------------------------------------------------%
        function uR = evaluateOnDomain(prop,model, state)
            [p, T] = model.getProps(state, 'pressure', 'Temperature');
            uR     = prop.evaluateFunctionOnDomainWithArguments(model.rock.uR, p, T); 
        end 
    end
    
end