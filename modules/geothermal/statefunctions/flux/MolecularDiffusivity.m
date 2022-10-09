classdef MolecularDiffusivity < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function md = MolecularDiffusivity(model, varargin)
            md@StateFunction(model, varargin{:});
            md = md.dependsOn({'Temperature', 'PoreVolume'}, 'PVTPropertyFunctions');
            md = md.dependsOn({'pressure'}, 'state');
            md.label = 'd_i';
        end
        
        %-----------------------------------------------------------------%
        function d = evaluateOnDomain(prop, model, state)
            % Get dependencies
            [p, T, pv] = model.getProps(state, 'pressure', 'Temperature', 'PoreVolume');
            % Compute tourtuosity and porosity
            tau  = prop.evaluateFunctionOnDomainWithArguments(model.rock.tau, p, T);
            poro = pv./model.G.cells.volumes;
            % Comoute diffusivity
            ncomp = model.getNumberOfComponents();
            d     = cell(ncomp,1);
            for i = 1:ncomp
                d{i} = model.Components{i}.molecularDiffusivity.*tau.*poro;
            end
        end
    end
    
end