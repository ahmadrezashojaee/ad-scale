classdef GeothermalFlowDiscretization < FlowDiscretization
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function props = GeothermalFlowDiscretization(model)
            % Inherit most of the state functions from FLuxDiscretization
            props = props@FlowDiscretization(model);
            % Set Thermal conductivity and transmissibility
            props = props.setStateFunction('RockThermalConductivity'  , RockThermalConductivity(model));
            props = props.setStateFunction('RockHeatTransmissibility' , DynamicTransmissibility(model, 'RockThermalConductivity'));
            props = props.setStateFunction('FluidThermalConductivity' , FluidThermalConductivity(model));
            props = props.setStateFunction('FluidHeatTransmissibility', DynamicTransmissibility(model, 'FluidThermalConductivity'));
            % Set conductive and advective heat fluxes
            props = props.setStateFunction('ConductiveHeatFlux', ConductiveHeatFlux(model));
            props = props.setStateFunction('AdvectiveHeatFlux' , AdvectiveHeatFlux(model));
            % Set molecular diffusion flux
            props = props.setStateFunction('MolecularDiffusivity'       , MolecularDiffusivity(model));
            props = props.setStateFunction('MolecularTransmissibility'  , DynamicTransmissibility(model, 'MolecularDiffusivity'));
            props = props.setStateFunction('ComponentTotalDiffusiveFlux', ComponentTotalDiffusiveFlux(model));
        end
        
        %-----------------------------------------------------------------%
        function [acc, flux, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            [acc, flux, names, types] = componentConservationEquations@FlowDiscretization(fd, model, state, state0, dt);
            % Conductive and advective heat flux
            flowState = fd.buildFlowState(model, state, state0, dt);
            diffFlux = model.getProp(flowState, 'ComponentTotalDiffusiveFlux');
            for i = 1:numel(flux)
                flux{i} = flux{i} + diffFlux{i};
            end
        end
        
        %-----------------------------------------------------------------%
        function [acc, flux, name, type] = energyConservationEquation(fd, model, state, state0, dt)
            nph = model.getNumberOfPhases;
            % Thermal energy in the fluid
            energyFluid  = model.getProps(state , 'PhaseThermalEnergy');
            energyFluid0 = model.getProps(state0, 'PhaseThermalEnergy');
            % Thermal energy in the fluid at previous timestep
            energyRock  = model.getProps(state , 'RockThermalEnergy');
            energyRock0 = model.getProps(state0, 'RockThermalEnergy');
            % Conductive and advective heat flux
            flowState = fd.buildFlowState(model, state, state0, dt);
            heatFluxCond = model.getProps(flowState, 'ConductiveHeatFlux');
            heatFluxAdv  = model.getProps(flowState, 'AdvectiveHeatFlux');
            if ~iscell(energyFluid0)
                % May have been reduced to a double array
                energyFluid0 = {energyFluid0};
            end
            % Add up accumulation
            acc = (energyRock  - energyRock0 )./dt;
            for i = 1:nph
                acc = acc + (energyFluid{i} - energyFluid0{i})./dt;
            end
            % Add up fluxes
            flux = heatFluxCond;
            for i = 1:nph
                flux = flux + heatFluxAdv{i};
            end
            % Convert to cell arrays to comply with AD framework
            acc  = {acc};
            flux = {flux};
            % Name and type
            name = 'energy';
            type = 'cell';
        end
    end
    
end