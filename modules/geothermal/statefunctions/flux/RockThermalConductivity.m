classdef RockThermalConductivity < StateFunction
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RockThermalConductivity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures', 'Temperature', 'PoreVolume'}, 'PVTPropertyFunctions');
            gp.label = '\Lambda_R';
        end
        
        %-----------------------------------------------------------------%
        function lambdaR = evaluateOnDomain(prop, model, state)
            [p, T, pv] = model.getProps(state, 'PhasePressures', 'Temperature', 'PoreVolume');
            lambdaR = prop.evaluateFunctionOnDomainWithArguments(model.rock.lambdaR, p, T);
            vol = model.G.cells.volumes;
            lambdaR = lambdaR.*(vol - pv)./vol;
        end
    end
    
end
