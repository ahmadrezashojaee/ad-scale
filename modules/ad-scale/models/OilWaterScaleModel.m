classdef OilWaterScaleModel < TwoPhaseOilWaterModel
    % Oil/water/polymer system
    %
    %
    % SYNOPSIS:
    %   model = OilWaterSclaeModel(G, rock, fluid, varargin)
    %
    % DESCRIPTION:
    %   Two phase model with scale deposition. A description of the scale model
    %   that is implemented here can be found in the directory ad-eor/docs .
    %
    % PARAMETERS:
    %   G        - Grid
    %   rock     - Rock structure
    %   fluid    - Fluid structure
    %   varargin - optional parameter
    %
    % RETURNS:
    %   class instance
    %
    % EXAMPLE:
    %
    % SEE ALSO: ThreePhaseBlackOilPolymerModel
    %

    properties
    
    Scale
    Sulfate
    Calcium
    Strontium
    Barium
    Sodium
    Magnesium
    Carbonate
    Chlorine
    iphreeqc
    prop
    initial_props
    pvCorrector
    
    end

    methods
        function model = OilWaterScaleModel(G, rock, fluid, varargin)

            model = model@TwoPhaseOilWaterModel(G, rock, fluid);

            % This is the model parameters for oil/water/scale
            model.oil = true;
            model.gas = false;
            model.water = true;
            model.Scale = true;
            model.Sulfate = true;
            model.Calcium = true;
            model.Strontium = true;
            model.Barium = true;
            model.Sodium = true;
            model.Magnesium = true;
            model.Carbonate = true;
            model.Chlorine = true;
            model.Components={SulfateComponent(), CalciumComponent(),...
                               StrontiumComponent(), BariumComponent(),...
                               SodiumComponent(), MagnesiumComponent(),...
                               CarbonateComponent(), ChlorineComponent()};
            model.iphreeqc = actxserver('IPhreeqcCOM.Object');
%             model.iphreeqc.LoadDatabase('C:\Program Files\USGS\IPhreeqcCOM 3.6.2-15100\database\phreeqc.dat');
            model.iphreeqc.LoadDatabase('phreeqc.dat'); 
            model.toleranceMB=1e-7;
            model.toleranceCNV=5e-3;
            model.FacilityModel.toleranceWellRate=1e-3;
            model.FacilityModel.toleranceWellMS=1e-3;
            model.nonlinearTolerance=1e-4;
            model = merge_options(model, varargin{:});
        end
        
        
        function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);
            % sulfate, calcium, strontium, barium must be present
            nc = model.G.cells.num;
            model.checkProperty(state, 'sulfate', [nc, 1], [1, 2]);
            model.checkProperty(state, 'calcium', [nc, 1], [1, 2]);
            model.checkProperty(state, 'strontium', [nc, 1], [1, 2]);
            model.checkProperty(state, 'barium', [nc, 1], [1, 2]);
            model.checkProperty(state, 'sodium', [nc, 1], [1, 2]);
            model.checkProperty(state, 'magnesium', [nc, 1], [1, 2]);
            model.checkProperty(state, 'carbonate', [nc, 1], [1, 2]);
            model.checkProperty(state, 'chlorine', [nc, 1], [1, 2]);
            
        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterScale(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            index = ':';
            switch(lower(name))
                case {'sulfate'}
                    fn = 'c_SO4';
                case {'calcium'}
                    fn = 'c_Ca';
                case {'strontium'}
                    fn = 'c_Sr';
                case {'barium'}
                    fn = 'c_Ba';
                case {'sodium'}
                    fn = 'c_Na';
                case {'magnesium'}
                    fn = 'c_Mg';
                case {'carbonate'}
                    fn = 'c_CO3';
                case {'chlorine'}
                    fn = 'c_Cl';
                case {'qw_so4'}
                    fn = 'qW_SO4';
                case {'qw_ca'}
                    fn = 'qW_Ca';
                case {'qw_sr'}
                    fn = 'qW_Sr';
                case {'qw_ba'}
                    fn = 'qW_Ba';
                case {'qw_na'}
                    fn = 'qW_Na';
                case {'qw_mg'}
                    fn = 'qW_Mg';
                case {'qw_co3'}
                    fn = 'qW_CO3';
                case {'qw_cl'}
                    fn = 'qW_Cl';
                otherwise
                    [fn, index] = getVariableField@ThreePhaseBlackOilModel(...
                                    model, name, varargin{:});
            end
        end

        function names = getComponentNames(model)
            names = getComponentNames@TwoPhaseOilWaterModel(model);
            if model.Sulfate
                names{end+1} = 'sulfate';
            end
            if model.Calcium
                names{end+1} = 'calcium';
            end
            if model.Strontium
                names{end+1} = 'strontium';
            end
            if model.Barium
                names{end+1} = 'barium';
            end
            if model.Sodium
                names{end+1} = 'sodium';
            end
            if model.Magnesium
                names{end+1} = 'magnesium';
            end
            if model.Carbonate
                names{end+1} = 'carbonate';
            end
            if model.Chlorine
                names{end+1} = 'chlorine';
            end
        end
   
        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        % For a given component conservation equation, compute and add in
        % source terms for a specific source/bc where the fluxes have
        % already been computed.
        %
        % PARAMETERS:
        %
        %   model  - (Base class, automatic)
        %
        %   cname  - Name of the component. Must be a property known to the
        %            model itself through `getProp` and `getVariableField`.
        %
        %   eq     - Equation where the source terms are to be added. Should
        %            be one value per cell in the simulation grid (model.G)
        %            so that the src.sourceCells is meaningful.
        %
        %   component - Cell-wise values of the component in question. Used
        %               for outflow source terms only.
        %
        %   src    - Source struct containing fields for fluxes etc. Should
        %            be constructed from force and the current reservoir
        %            state by `computeSourcesAndBoundaryConditionsAD`.
        %
        %   force  - Force struct used to produce src. Should contain the
        %            field defining the component in question, so that the
        %            inflow of the component through the boundary condition
        %            or source terms can accurately by estimated.
            if isempty(force)
                return
            end
            c = model.getProp(force, cname);
            cells = src.sourceCells;
            switch lower(cname)
              case {'sulfate', 'calcium', 'strontium', 'barium', 'sodium', 'magnesium', 'carbonate', 'chlorine'}
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
              otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
            end
            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;
        end
        
        function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@ThreePhaseBlackOilModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.Sulfate || model.Calcium || model.strontium || model.Barium
                assert(model.water, 'Scale deposition requires a water phase.');
                f = model.fluid;

                % Water is always first
                wix = 1;
                %cqWs is water injection rate m^3/s
                cqWs = qMass{wix}./f.rhoWS; % connection volume flux at surface condition
            end

            if model.Sulfate
                if well.isInjector()
                    concWells = model.getProp(well.W, 'sulfate');
                else
                    pixs = strcmpi(model.getComponentNames(), 'sulfate');
                    concWells = packed.components{pixs};
                end
                cq_SO4 = concWells.*cqWs; %sulfate rate: mass rate = cq_SO4*f.rhoWS
                qw_SO4 = packed.extravars{strcmpi(packed.extravars_names, 'qw_SO4')};

                compEqs{end+1} = qw_SO4 - sum(concWells.*cqWs);
                compSrc{end+1} = cq_SO4; %sulfate source
                eqNames{end+1} = 'sulfateWells';
            end
            if model.Calcium
                if well.isInjector()
                    concWells = model.getProp(well.W, 'calcium');
                else
                    pixs = strcmpi(model.getComponentNames(), 'calcium');
                    concWells = packed.components{pixs};
                end
                cq_Ca = concWells.*cqWs;
                qw_Ca = packed.extravars{strcmpi(packed.extravars_names, 'qw_Ca')};

                compEqs{end+1} = qw_Ca - sum(concWells.*cqWs);
                compSrc{end+1} = cq_Ca;
                eqNames{end+1} = 'calciumWells';
            end
            if model.Strontium
                if well.isInjector()
                    concWells = model.getProp(well.W, 'strontium');
                else
                    pixs = strcmpi(model.getComponentNames(), 'strontium');
                    concWells = packed.components{pixs};
                end
                cq_Sr = concWells.*cqWs;
                qw_Sr = packed.extravars{strcmpi(packed.extravars_names, 'qw_Sr')};

                compEqs{end+1} = qw_Sr - sum(concWells.*cqWs);
                compSrc{end+1} = cq_Sr;
                eqNames{end+1} = 'strontiumWells';
            end

            if model.Barium
                if well.isInjector()
                    concWells = model.getProp(well.W, 'barium');
                else
                    pixs = strcmpi(model.getComponentNames(), 'barium');
                    concWells = packed.components{pixs};
                end
                cq_Ba = concWells.*cqWs;
                qw_Ba = packed.extravars{strcmpi(packed.extravars_names, 'qw_Ba')};

                compEqs{end+1} = qw_Ba - sum(concWells.*cqWs);
                compSrc{end+1} = cq_Ba;
                eqNames{end+1} = 'bariumWells';
            end
            
            if model.Sodium
                if well.isInjector()
                    concWells = model.getProp(well.W, 'sodium');
                else
                    pixs = strcmpi(model.getComponentNames(), 'sodium');
                    concWells = packed.components{pixs};
                end
                cq_Na = concWells.*cqWs;
                qw_Na = packed.extravars{strcmpi(packed.extravars_names, 'qw_Na')};

                compEqs{end+1} = qw_Na - sum(concWells.*cqWs);
                compSrc{end+1} = cq_Na;
                eqNames{end+1} = 'sodiumWells';
            end
            
            if model.Magnesium
                if well.isInjector()
                    concWells = model.getProp(well.W, 'magnesium');
                else
                    pixs = strcmpi(model.getComponentNames(), 'magnesium');
                    concWells = packed.components{pixs};
                end
                cq_Mg = concWells.*cqWs;
                qw_Mg = packed.extravars{strcmpi(packed.extravars_names, 'qw_Mg')};

                compEqs{end+1} = qw_Mg - sum(concWells.*cqWs);
                compSrc{end+1} = cq_Mg;
                eqNames{end+1} = 'magnesiumWells';
            end
            
            if model.Carbonate
                if well.isInjector()
                    concWells = model.getProp(well.W, 'carbonate');
                else
                    pixs = strcmpi(model.getComponentNames(), 'carbonate');
                    concWells = packed.components{pixs};
                end
                cq_CO3 = concWells.*cqWs;
                qw_CO3 = packed.extravars{strcmpi(packed.extravars_names, 'qw_CO3')};

                compEqs{end+1} = qw_CO3 - sum(concWells.*cqWs);
                compSrc{end+1} = cq_CO3;
                eqNames{end+1} = 'carbonateWells';
            end
            
            if model.Chlorine
                if well.isInjector()
                    concWells = model.getProp(well.W, 'chlorine');
                else
                    pixs = strcmpi(model.getComponentNames(), 'chlorine');
                    concWells = packed.components{pixs};
                end
                cq_Cl = concWells.*cqWs;
                qw_Cl = packed.extravars{strcmpi(packed.extravars_names, 'qw_Cl')};

                compEqs{end+1} = qw_Cl - sum(concWells.*cqWs);
                compSrc{end+1} = cq_Cl;
                eqNames{end+1} = 'chlorineWells';
            end
        end
        
        function names = getExtraWellPrimaryVariableNames(model)
            names = getExtraWellPrimaryVariableNames@ThreePhaseBlackOilModel(model);
            if model.Sulfate
                names{end+1} = 'qW_SO4';
            end
            if model.Calcium
                names{end+1} = 'qW_Ca';
            end
            if model.Strontium
                names{end+1} = 'qW_Sr';
            end
            if model.Barium
                names{end+1} = 'qW_Ba';
            end
            if model.Sodium
                names{end+1} = 'qW_Na';
            end
            if model.Magnesium
                names{end+1} = 'qW_Mg';
            end
            if model.Carbonate
                names{end+1} = 'qW_CO3';
            end
            if model.Chlorine
                names{end+1} = 'qW_Cl';
            end
        end
        
        function [names, types] = getExtraWellEquationNames(model)
            [names, types] = getExtraWellEquationNames@ThreePhaseBlackOilModel(model);
            if model.Sulfate
                names{end+1} = 'sulfateWells';
                types{end+1} = 'perf';
            end
            if model.Calcium
                names{end+1} = 'calciumWells';
                types{end+1} = 'perf';
            end
            if model.Strontium
                names{end+1} = 'strontiumWells';
                types{end+1} = 'perf';
            end
            if model.Barium
                names{end+1} = 'bariumWells';
                types{end+1} = 'perf';
            end
            if model.Sodium
                names{end+1} = 'sodiumWells';
                types{end+1} = 'perf';
            end
            if model.Magnesium
                names{end+1} = 'magnesiumWells';
                types{end+1} = 'perf';
            end
            if model.Carbonate
                names{end+1} = 'carbonateWells';
                types{end+1} = 'perf';
            end
            if model.Chlorine
                names{end+1} = 'chlorineWells';
                types{end+1} = 'perf';
            end
        end
        
    end
end


%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
