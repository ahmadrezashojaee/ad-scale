classdef Enthalpy < StateFunction
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = Enthalpy(model, varargin)
            gp@StateFunction(model, varargin{:});
            switch model.thermalFormulation
                case 'enthalpy'
                    % Enthalpy is the primary variable
                case 'temperature'
                    % Temperature is the primary variable
                    gp = gp.dependsOn({'PhasePressures', 'Density', 'PhaseInternalEnergy'});
            end
            gp.label = 'h_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function h = evaluateOnDomain(prop, model, state)
            switch model.thermalFormulation
                case 'enthalpy'
                    error('Not implemented yet');
                case 'temperature'
                    [p, rho, u] = prop.getEvaluatedDependencies(state, 'PhasePressures'     , ...
                                                                       'Density'            , ...
                                                                       'PhaseInternalEnergy');
                    nph = model.getNumberOfPhases();
                    h = cell(1, nph);
                    for i = 1:nph
                        h{i} = u{i} + p{i}./rho{i};
                    end
            end
        end
    end

end