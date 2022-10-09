classdef PhaseInternalEnergy < StateFunction
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = PhaseInternalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures', 'Temperature', 'ComponentPhaseMassFractions'});
            gp.label = 'u_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function u = evaluateOnDomain(prop,model, state)
            [p, T, X] = prop.getEvaluatedDependencies(state, 'PhasePressures'             , ...
                                                             'Temperature'                , ...
                                                             'ComponentPhaseMassFractions');
            cnames = model.getComponentNames();
            ix = strcmpi(cnames, 'NaCl');
            if any(ix), X = X{ix}; else, X = []; end
            nph = model.getNumberOfPhases();
            u   = cell(1, nph);
            if model.water
                wix    = model.getPhaseIndex('W');
                pw     = p{wix};
                u{wix} = prop.evaluateFluid(model, 'uW', pw, T, X);
            end
            if model.oil
                oix    = model.getPhaseIndex('O');
                po     = p{oix};
                u{oix} = prop.evaluateFluid('uO', po, T);
            end
            if model.gas
                gix    = model.getPhaseIndex('G');
                pg     = p{gix};
                u{gix} = prop.evaluateFluid('uG', pg, T);
            end
        end       
    end
end

