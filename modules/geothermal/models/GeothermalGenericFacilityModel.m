classdef GeothermalGenericFacilityModel < GenericFacilityModel
   
    properties
        thermalFormulation = 'temperature';
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, useDefaults)
            if nargin < 2
                useDefaults = isempty(model.FacilityFlowDiscretization);
            end
            model = setupStateFunctionGroupings@GenericFacilityModel(model, useDefaults);
            ffd = model.FacilityFlowDiscretization;
            ffd = ffd.setStateFunction('AdvectiveHeatFlux' , WellAdvectiveHeatFlux(model) );
            ffd = ffd.setStateFunction('ConductiveHeatFlux', WellConductiveHeatFlux(model));
            model.FacilityFlowDiscretization = ffd;
        end
        
        %-----------------------------------------------------------------%
        function names = getBasicPrimaryVariableNames(model)
            names = getBasicPrimaryVariableNames@GenericFacilityModel(model);
            if strcmpi(model.thermalFormulation, 'none')
                return
            end
            if model.ReservoirModel.thermal
                names = [names, 'well_temperature'];
            end
        end
        
        %-----------------------------------------------------------------%
        function [variables, names, map] = getBasicPrimaryVariables(model, wellSol)
            [variables, names, map] = getBasicPrimaryVariables@GenericFacilityModel(model, wellSol);
            if model.ReservoirModel.thermal
                isBHP  = strcmpi(names, 'bhp');
                isTemp = strcmpi(names, 'well_temperature');
                isRate = ~isBHP & ~isTemp;
                map = struct('isBHP', isBHP, 'isRate', isRate, 'isTemp', isTemp);
            end
        end
        
        %-----------------------------------------------------------------%
        function src = getEnergySources(facility, state)
            map = facility.getProps(state, 'FacilityWellMapping');
            if isempty(map.W)
                val = [];
            else
                val = facility.getProps(state, 'ConductiveHeatFlux');
                qa = facility.getProps(state, 'AdvectiveHeatFlux');
                for i = 1:numel(qa)
                    val = val + qa{i};
                end
            end
            src = struct('value', {{val}}, 'cells', map.cells);
        end
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            [eqs, names, types, state] = getModelEquations@GenericFacilityModel(model, state0, state, dt, drivingForces);
            isTemp = strcmpi(state.FacilityState.names, 'well_temperature');
            if strcmpi(model.thermalFormulation, 'none') || ~any(isTemp)
                return
            end
            T  = model.getWellTemperature(state);
            Tw = state.FacilityState.primaryVariables{isTemp};
            eqs = [eqs, {Tw - T}];
            names = [names, {'temperatureWells'}];
            types = [types, {'perf'}];
        end
        
        %-----------------------------------------------------------------%
        function T = getWellTemperature(model, state)
            [map, q] = model.getProps(state, 'FacilityWellMapping', 'PhaseFlux');
            [T, pv] = model.ReservoirModel.getProps(state, 'Temperature', 'PoreVolume');
            T  = T(map.cells);
            pv = pv(map.cells);
            qT = sum(value(q),2);
            injector = qT > 0;
            if any(injector)
                TWell = vertcat(map.W.T);
                TWell = TWell(map.perf2well);
                T(injector) = TWell(injector);
            end
            pvt = map.perforationSum*pv;
            T = (map.perforationSum*(pv.*T))./pvt;
        end

        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@GenericFacilityModel(model, state0, state, dt, drivingForces);
            map = model.getProp(state, 'FacilityWellMapping');
            if strcmpi(model.thermalFormulation, 'none')
                Tw = num2cell(model.getWellTemperature(state));
                [state.wellSol(map.active).T] = deal(Tw{:});
            end
            if numel(map.active) < numel(state.wellSol)
                inactive = true(numel(state.wellSol),1);
                inactive(map.active) = false;
                if isfield(state0, 'wellSol')
                    [state.wellSol(inactive).T] = deal(state0.wellSol(inactive).T);
                end
            end
        end
    end
    
end