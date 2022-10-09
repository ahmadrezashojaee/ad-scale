classdef ThermalDensity < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = ThermalDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures', 'Temperature', 'ComponentPhaseMassFractions'});
            gp.label = '\rho_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function rho = evaluateOnDomain(prop, model, state)
            [p, T, X] = prop.getEvaluatedDependencies(state, 'PhasePressures'             , ...
                                                             'Temperature'                , ...
                                                             'ComponentPhaseMassFractions');
            cnames = model.getComponentNames();
            ix = strcmpi(cnames, 'NaCl');
            if any(ix), X = X{ix}; else, X = []; end
            phases = model.getPhaseNames();
            rho = cell(1, numel(phases));
            for i = 1:numel(phases)
                rho{i} = prop.evaluateFluid(model, ['rho', phases(i)], p{i}, T, X);
            end
        end
    end
    
end

