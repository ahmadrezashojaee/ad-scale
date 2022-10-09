classdef GeothermalModel < ReservoirModel & GenericReservoirModel
   
    properties
        % Thermal properties
        thermal            = true;
        thermalFormulation = 'temperature';
        % Component properties
        compFluid       % Component data, similar to cmopositional fluid
        dxMaxAbs = 0.1; % Maximum allowed update to molar fractions
        % Maximum and minimum allowed temperature
        minimumTemperature = -inf;
        maximumTemperature = inf;
        dTMaxRel = 0.2;
        dTMaxAbs = Inf;
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = GeothermalModel(G, rock, fluid, compFluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            model.water = true;
            if nargin < 4 || isempty(compFluid)
                % Default is a one-component system with H20 only
                compFluid = CompositionalBrineFluid({'H2O'}           , ...
                                                    18.015281*gram/mol, ...
                                                    0                 );
            end
            model.compFluid = compFluid;
            model = merge_options(model, varargin{:});
            model.OutputStateFunctions = {'ComponentTotalMass', 'PhaseThermalEnergy', 'RockThermalEnergy'};
            assert(strcmpi(model.thermalFormulation, 'temperature'), ...
                'GeothermalModel currently only supports temperature formulation');
        end
        
        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            if isempty(model.FacilityModel) || ~isa(model.FacilityModel, 'GeothermalGenericFacilityModel')
                model.FacilityModel = GeothermalGenericFacilityModel(model);
            end
            if isempty(model.Components)
                names = model.getComponentNames();
                nc = numel(names);
                for i = 1:nc
                    name        = names{i};
                    molarMass   = model.compFluid.molarMass(i);
                    diffusivity = model.compFluid.molecularDiffusivity(i);
                    c = BrineComponent(name, molarMass, diffusivity, i);
                    model.Components{i} = c;
                end
            end
            model = validateModel@ReservoirModel(model, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function ctrl = validateDrivingForces(model, ctrl, index)
            ctrl = validateDrivingForces@ReservoirModel(model, ctrl, index);
            if isfield(ctrl, 'bc') && ~isempty(ctrl.bc)
                % Check that bcs have thermal fields
                assert((isfield(ctrl.bc, 'T'    ) && numel(ctrl.bc.T)     == numel(ctrl.bc.face)) || ...
                       (isfield(ctrl.bc, 'Hflux') && numel(ctrl.bc.Hflux) == numel(ctrl.bc.Hflux)),  ...
                       'bc must have a given temperature (T) or heat flux (Hflux)'                );
            end
            if isfield(ctrl, 'W') && ~isempty(ctrl.W)
                % Check that wells have a prescribed temperature
                assert(all(arrayfun(@(w) isfield(w, 'T'), ctrl.W))  , ...
                       'All wells must have a temperature field (T)');
            end
        end
        
        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@ReservoirModel(model, varargin{:});
            % PVT properties
            pvt = model.PVTPropertyFunctions;
            pvt = pvt.setStateFunction('ShrinkageFactors'           , DensityDerivedShrinkageFactors(model));
            pvt = pvt.setStateFunction('Density'                    , ThermalDensity(model));
            pvt = pvt.setStateFunction('ComponentPhaseMassFractions', ComponentPhaseMassFractionsBrine(model));
            pvt = pvt.setStateFunction('ComponentPhaseMoleFractions', ComponentPhaseMoleFractionsBrine(model));
            pvt = pvt.setStateFunction('Viscosity'                  , ThermalViscosity(model));
            % Temperature and enthalpy
            pvt = pvt.setStateFunction('Temperature'        , Temperature(model));
            pvt = pvt.setStateFunction('Enthalpy'           , Enthalpy(model));
            pvt = pvt.setStateFunction('PhaseInternalEnergy', PhaseInternalEnergy(model));
            pvt = pvt.setStateFunction('RockInternalEnergy' , RockInternalEnergy(model));
            % Thermal energy in fluid phases
            pvt = pvt.setStateFunction('PhaseMass'         , PhaseMass(model));
            pvt = pvt.setStateFunction('PhaseThermalEnergy', PhaseThermalEnergy(model));
            % Thermal energy in rock
            pvt = pvt.setStateFunction('RockMass'         , RockMass(model));
            pvt = pvt.setStateFunction('RockDensity'      , RockDensity(model));
            pvt = pvt.setStateFunction('RockThermalEnergy', RockThermalEnergy(model));
            % Replace
            model.PVTPropertyFunctions = pvt;
            % Make flux discretization
            warning('Assuming default flux discretization');
            fd = GeothermalFlowDiscretization(model);
            model.FlowDiscretization = fd; 
        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            state = validateState@ReservoirModel(model, state);
            if ~isfield(state, 'T')
                T0 = (273.15 + 20)*Kelvin;
                state.T = repmat(T0, model.G.cells.num, 1);
            end
            if ~model.thermal
                return
            end
            if isfield(state, 'enthalphy')
                warning('Enthalpy given in initial state. I''ll assume this is consistent with the EOS and carry on');
                return
            end
            if strcmpi(model.thermalFormulation, 'enthalpy')
                [p, T, rho] = model.getProps(state, 'PhasePressures', 'temperature', 'Density');
                state.enthalpy = T.*model.fluid.cpW(p{1}, T) + p{1}./rho{1};
            end
            ncomp = model.getNumberOfComponents();
            if isfield(state, 'components')
                assert(all(size(state.components) == [model.G.cells.num, ncomp]));
            else
                state.components = zeros(model.G.cells.num, ncomp);
                state.components(:,1) = 1;
            end
        end
        
        %-----------------------------------------------------------------%
        function [vararg, control] = getDrivingForces(model, control)
            [vararg, control] = getDrivingForces@ReservoirModel(model, control);
            ix = find(strcmpi(vararg, 'bc'));
            if ~isempty(ix) && ~isempty(vararg{ix+1})
                bc = vararg{ix+1};
                if ~isfield(bc, 'components')
                    assert(numel(model.compFluid.names) == 1, ...
                        ['Model has more than one component - please ', ...
                         'provide component field to bc (bc.cmomponents)']);
                    bc.components = ones(numel(bc.face),1);
                    vararg{ix+1}  = bc;
                end
            end
            ix = find(strcmpi(vararg, 'W'));
            if ~isempty(ix) && ~isempty(vararg{ix+1})
                W = vararg{ix+1};
                if isempty(W(1).components)
                    assert(numel(model.compFluid.names) == 1, ...
                        ['Model has more than one component - ...    ', ...
                         'please provide component field to well (W.components)']);
                    [W.components] = deal(1);
                    vararg{ix+1} = W;
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function [vars, names, origin] = getPrimaryVariables(model, state)
            [p, x] = model.getProps(state, 'pressure', 'components');
            x = ensureMinimumFraction(x, 1e-8);
            x = expandMatrixToCell(x);
            cnames = model.getComponentNames();
            vars   = [p, x(:, 2:end)];
            names  = [{'pressure'}, cnames(2:end)];
            if model.thermal
                thName = model.thermalFormulation;
                thVar  = model.getProps(state, thName);
                vars   = [vars , thVar ];
                names  = [names, thName];
            end
            origin = repmat({class(model)}, 1, numel(cnames) + model.thermal);
            if ~isempty(model.FacilityModel)
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
                vars = [vars, v];
                names = [names, n];
                origin = [origin, o];
            end
        end
        
        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)
            isP = strcmp(names, 'pressure');
            state = model.setProp(state, 'pressure', vars{isP});
            removed = isP;
            
            cnames = model.getComponentNames();
            ncomp = numel(cnames);
            x = cell(1, ncomp);
            x_end = ones(model.G.cells.num, 1);
            for i = 1:ncomp
                name = cnames{i};
                sub = strcmp(names, name);
                if any(sub)
                    x{i} = vars{sub};
                    x_end = x_end - x{i};
                    removed(sub) = true;
                else
                    fill = i;
                end
            end
            x{fill} = x_end;
            state = model.setProp(state, 'components', x);
            if ~isempty(model.FacilityModel)
                % Select facility model variables and pass them off to attached
                % class.
                fm = class(model.FacilityModel);
                isF = strcmp(origin, fm);
                state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                removed = removed | isF;
            end
            % Set up state with remaining variables
            state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));
        end

        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            [eqs, flux, names, types] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Add sources
            eqs = model.insertSources(eqs, src);
            % Assemble equations
            for i = 1:numel(eqs)
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end
            if model.thermal
                [eeqs, eflux, enames, etypes] = model.FlowDiscretization.energyConservationEquation(model, state, state0, dt);
                % Add sources
                esrc = model.FacilityModel.getEnergySources(state);
                eeqs = model.insertSources(eeqs, esrc);
                % Assemble equations
                eeqs{1} = model.operators.AccDiv(eeqs{1}, eflux{1});
            else
                [eeqs, enames, etypes] = deal([]);
            end
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            % Concatenate
            eqs   = [eqs  , eeqs  , weqs  ];
            names = [names, enames, wnames];
            types = [types, etypes, wtypes];
            % Add in boundary conditions
            eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, drivingForces);
            % Scaling
            scale = model.getEquationScaling(eqs, names, state0, dt);
            for i = 1:numel(eqs)
                if ~isempty(scale{i})
                    eqs{i} = eqs{i}.*scale{i};
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function scale = getEquationScaling(model, eqs, names, state0, dt)
            scale = cell(numel(eqs),1);
            cnames = model.getComponentNames();
            [cmass, energyRock, energyPhase] = model.getProps(state0, 'ComponentTotalMass', ...
                                                                     'RockThermalEnergy' , ...
                                                                     'PhaseThermalEnergy');
            cmass = value(cmass); energyRock = value(energyRock); energyPhase = value(energyPhase);
            if ~iscell(cmass), cmass = {cmass}; end
            ncomp = model.getNumberOfPhases();
            mass = 0;
            for i = 1:ncomp
                mass = mass + cmass{i};
            end
            scaleMass = dt./mass;
            for n = cnames
               ix = strcmpi(n{1}, names);
               if ~any(ix)
                   continue
               end
               scale{ix} = scaleMass;
            end
            ix = strcmpi(names, 'energy');
            if ~isempty(ix)
                scaleEnergy = dt./(energyRock + sum(energyPhase, 2));
                scale{ix} = scaleEnergy;
            end
        end
        
        %-----------------------------------------------------------------%
        function massFraction = getMassFraction(model, moleFraction)
            % Convert molar fraction to mass fraction
            if iscell(moleFraction)
                ncomp = numel(moleFraction);
                mass = cell(1, ncomp);
                totMass = 0;
                for i = 1:ncomp
                    mi = moleFraction{i};
                    if ~isempty(mi)
                        mass{i} = model.Components{i}.molarMass.*moleFraction{i};
                        totMass = totMass + mass{i};
                    end
                end
                massFraction = cell(size(mass));
                for i = 1:ncomp
                    if ~isempty(mass{i})
                        massFraction{i} = mass{i}./totMass;
                    end
                end
            else
                molarMass = reshape(cellfun(@(c) c.molarMass, model.Components), 1, []);
                mass = bsxfun(@times, moleFraction, molarMass);
                massFraction = bsxfun(@rdivide, mass, sum(mass, 2));
            end
        end
        
        %-----------------------------------------------------------------%
        function moleFraction = getMoleFraction(model, massFraction)
            % Convert molar fraction to mass fraction
            if iscell(massFraction)
                assert(~isa(massFraction{1}, 'ADI'));
                massFraction = horzcat(massFraction{:});
            end
            model = model.validateModel();
            mass = massFraction;
            molarMass = reshape(cellfun(@(c) c.molarMass, model.Components), 1, []);
            mole = bsxfun(@rdivide, mass, molarMass);
            moleFraction = bsxfun(@rdivide, mole, sum(mole, 2));
        end
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            nm = lower(name);
            switch nm
                % Make sure we go through state functions instead of state
                % for the derived variable
                case {'enthalpy', 'h'}
                    if ~strcmpi(model.thermalFormulation, 'enthalpy')
                        nm = '';
                    end
                case {'temperature', 't'}
                    if ~strcmpi(model.thermalFormulation, 'temperature')
                        nm = '';
                    end
            end
            switch nm
                case {'x', 'molefractions', 'components'}
                    % Liquid phase mole fraction
                    fn = 'components';
                    index = ':';
                case {'enthalpy', 'h'}
                    fn    = 'enthalpy';
                    index = ':';
                case {'temperature', 't'}
                    fn    = 'T';
                    index = ':';
                otherwise
                    cnames = model.getComponentNames();
                    sub = strcmpi(cnames, name);
                    if any(sub)
                        fn    = 'components';
                        index = find(sub);
                    else
                        % This will throw an error for us
                        [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
                    end
            end
        end
        
        %-----------------------------------------------------------------%
        function names = getComponentNames(model)
            names  = getComponentNames@ReservoirModel(model);
            names  = horzcat(names,  model.compFluid.names);
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dz, drivingForces)
            vars0 = problem.primaryVariables;
            vars = vars0;
            % Temperature
            state = model.updateStateFromIncrement(state, dz, problem, 'temperature', model.dTMaxRel, model.dTMaxAbs);
            state = model.capProperty(state, 'temperature', model.minimumTemperature, model.maximumTemperature);
            [vars, removed] = model.stripVars(vars, 'temperature');
            % Components
            cnames = model.getComponentNames();
            ncomp = numel(cnames);
            ok = false(ncomp, 1);
            x = state.components;
            rm = 0;
            for i = 1:ncomp
                name = lower(cnames{i});
                cix = strcmpi(vars0, name);
                if any(cix)
                    x0 = x(:, i);
                    dx = dz{cix};
                    if isfinite(model.dxMaxAbs)
                        dx = sign(dx).*min(abs(dx), model.dxMaxAbs);
                    end
                    x(:, i) = min(max(x0 + dx, 0), 1);
                    ok(i) = true;
                    [vars, ix] = model.stripVars(vars, {name});
                    removed(~removed) = removed(~removed) | ix;
                    rm = rm - (x(:, i) - x0);
                end
            end
            if any(ok)
                % We had components as active variables somehow
                assert(nnz(~ok) == 1)
                x(:, ~ok) = min(max(x(:, ~ok) + rm, 0), 1);
                x = bsxfun(@rdivide, x, sum(x, 2));
                state.components = x;
            else
                state.dx = zeros(1, ncomp);
            end

            % Parent class handles almost everything for us
            problem.primaryVariables = vars;
            dz(removed) = [];
            [state, report] = updateState@ReservoirModel(model, state, problem, dz, drivingForces);
            
            if problem.iterationNo == 1
                state.switched = false(model.G.cells.num, 1);
                state.switchCount = zeros(model.G.cells.num, 1);
            end
            state.components = ensureMinimumFraction(state.components, 1e-8);
        end
        
        %-----------------------------------------------------------------%
        function [eqs, state, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, forces)
             % Assemble equations and add in sources
            [pressures, sat, mob, rho, X] = model.getProps(state, 'PhasePressures'             , ...
                                                                  's'                          , ...
                                                                  'Mobility'                   , ...
                                                                  'Density'                    , ...
                                                                  'ComponentPhaseMassFractions');
            comps = cellfun(@(x) {x}, X, 'UniformOutput', false);
            [eqs, state, src] = addBoundaryConditionsAndSources@ReservoirModel(model, eqs, names, types, state, ...
                                                                      pressures, {sat}, mob, rho, {}, comps, forces);
            if ~model.thermal
                return
            end
            eix = strcmpi(names, 'energy');
            if ~isempty(src.bc.sourceCells)
                src = getHeatFluxBoundary(model, state, src, forces);
                q = src.bc.mapping*src.bc.heatFlux;
                eqs{eix}(src.bc.sourceCells) = eqs{eix}(src.bc.sourceCells) - q;
            end
        end
        
        %-----------------------------------------------------------------%
        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
            if isempty(force)
                return
            end
            cnames = model.getComponentNames();
            sub = strcmpi(cnames, cname);
            if any(sub)
                cells = src.sourceCells;
                x_bc = model.getProp(force, 'components');
                mf_bc = model.getMassFraction(x_bc);    
                massFractions = {mf_bc};
                qC = zeros(size(cells));
                nph = model.getNumberOfPhases;
                for ph = 1:nph
                    q_ph = src.phaseMass{ph};
                    inj = q_ph > 0;
                    qC = qC + ~inj.*component{ph}(cells).*q_ph ...
                            +  inj.*massFractions{ph}(:, sub).*q_ph;
                end
                if ~isempty(src.mapping)
                    qC = src.mapping*qC;
                end
                eq(cells) = eq(cells) - qC;
                src.components{end+1} = qC;
            else
                [eq, src] = addComponentContributions@ReservoirModel(model, cname, eq, component, src, force);
            end
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
           [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
           if model.extraStateOutput
               rho = model.getProps(state, 'Density');
               state.rho = horzcat(rho{:});
           end
           if model.outputFluxes
               [heatFluxAdv, heatFluxCond] = model.getProps(state, 'AdvectiveHeatFlux', 'ConductiveHeatFlux');
               state.heatFluxAdv  = horzcat(heatFluxAdv{:});
               state.heatFluxCond = heatFluxCond;
           end
        end
        
        % ----------------------------------------------------------------%
        function [model, state] = updateForChangedControls(model, state, forces)
            state0 = state;
            [model, state] = updateForChangedControls@ReservoirModel(model, state, forces);
            if isfield(state0, 'wellSol')
                inactive = ~vertcat(state.wellSol.status);
                [state.wellSol(inactive).T] = deal(state0.wellSol(inactive).T);
            end
        end
    end
    
end