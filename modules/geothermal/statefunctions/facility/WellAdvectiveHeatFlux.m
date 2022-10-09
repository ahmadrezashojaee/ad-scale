classdef WellAdvectiveHeatFlux < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = WellAdvectiveHeatFlux(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseFlux', 'FacilityWellMapping'});
            gp.label = 'q_{avd}';
        end
        
        %-----------------------------------------------------------------%
        function q = evaluateOnDomain(prop, model, state)
            [v, map] = prop.getEvaluatedDependencies(state, 'PhaseFlux', 'FacilityWellMapping');
            [h, rho] = model.ReservoirModel.getProps(state, 'Enthalpy', 'Density');
            h   = cellfun(@(h) h(map.cells), h, 'UniformOutput', false);
            rho = cellfun(@(rho) rho(map.cells), rho, 'UniformOutput', false);
            injector = map.isInjector(map.perf2well);            
            if any(injector)
                isbhp = strcmpi(state.FacilityState.names, 'bhp');
                if any(isbhp)
                    bhp = state.FacilityState.primaryVariables{isbhp}(map.perf2well);
                else
                    bhp = vertcat(wellSol(map.active).bhp);
                    bhp = bhp(map.perf2well);
                end
                T = [map.W.T]';
                T = T(map.perf2well);
                phases = model.ReservoirModel.getPhaseNames();
                for i = 1:numel(phases)
                    hWell = feval(model.ReservoirModel.fluid.(['h', phases(i)]), bhp, T);
                    h{i}(injector) = hWell(injector);
                end
            end
            q = cellfun(@(rho,v,h) rho.*v.*h, rho, v, h, 'UniformOutput', false);
        end 
    end
    
end