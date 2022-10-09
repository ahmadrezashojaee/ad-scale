classdef ThermalViscosity < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = ThermalViscosity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures'              , ...
                               'Temperature'                 , ...
                               'ComponentPhaseMassFractions'});
            gp.label = '\mu_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function mu = evaluateOnDomain(prop, model, state)
            [p, T, X] = prop.getEvaluatedDependencies(state, 'PhasePressures'             , ...
                                                             'Temperature'                , ...
                                                             'ComponentPhaseMassFractions');
            cnames = model.getComponentNames();
            ix = strcmpi(cnames, 'NaCl');
            if any(ix), X = X{ix}; else, X = []; end
            phases = model.getPhaseNames();
            mu = cell(1, numel(phases));
            for i = 1:numel(phases)
                mu{i} = prop.evaluateFluid(model, ['mu', phases(i)], p{i}, T, X);
            end
        end
    end
    
end

