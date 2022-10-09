function [problem, state] = equationsOilWaterScale(state0, state, model, dt, ...
                                                     drivingForces, varargin)
%
% SYNOPSIS:
%   function [problem, state] = equationsOilWaterScale(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: 
%   Assemble the linearized equations for an oil-water-Scale
%   system, computing both the residuals and the Jacobians. Returns the result as
%   an instance of the class LinearizedProblem which can be solved using instances
%   of LinearSolverAD.
%
% A description of the modeling equations can be found in the directory
% ad-eor/docs.
%
%
% PARAMETERS:
%   state0        - State at previous times-step
%   state         - State at current time-step
%   model         - Model instance
%   dt            - time-step
%   drivingForces - Driving forces (boundary conditions, wells, ...)
%   varargin      - optional parameters
%
% RETURNS:
%   problem - Instance of LinearizedProblem
%   state   - Updated state variable (fluxes, mobilities and more can be
%             stored, the wellSol structure is also updated in case of control switching)
%
% EXAMPLE:
%
% SEE ALSO: LinearizedProblem, LinearSolverAD, equationsOilWater, OilWaterPolymerModel
%

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


% Get linearized problem for oil/water/polymer system with black oil-style
% properties
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
op = model.operators;
fluid = model.fluid;

% Properties at current timestep
[p, sW, c_SO4,c_Ca,c_Sr,c_Ba, c_Na, c_Mg, c_CO3, c_Cl, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'sulfate', 'calcium','strontium','barium', 'sodium', 'magnesium', 'carbonate','chlorine', 'wellSol');

% Properties at previous timestep
[p0, sW0, c0_SO4,c0_Ca,c0_Sr,c0_Ba, c0_Na, c0_Mg, c0_CO3, c0_Cl, wellSol0] = model.getProps(state, 'pressure', 'water', ...
    'sulfate', 'calcium','strontium','barium', 'sodium', 'magnesium', 'carbonate','chlorine', 'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, c_SO4,c_Ca,c_Sr,c_Ba, c_Na, c_Mg, c_CO3, c_Cl, wellVars{:}] = ...
            model.AutoDiffBackend.initVariablesAD(p, sW, c_SO4,c_Ca,c_Sr,c_Ba, c_Na, c_Mg, c_CO3, c_Cl, wellVars{:});
        primaryVars = {'pressure', 'sW','sulfate', 'calcium','strontium','barium', 'sodium', 'magnesium', 'carbonate','chlorine', wellVarNames{:}};
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, sW0, c0_SO4,c0_Ca,c0_Sr,c0_Ba, c0_Na, c0_Mg, c0_CO3, c0_Cl, wellVars0{:}] = ...
            model.AutoDiffBackend.initVariablesAD(p0, sW0, c0_SO4, c0_Ca, c0_Sr, c0_Ba, c0_Na, c0_Mg, c0_CO3, c0_Cl, wellVars0{:}); %#ok
        primaryVars = {'pressure', 'sW','sulfate', 'calcium','strontium','barium', 'sodium', 'magnesium', 'carbonate','chlorine'};
    end
else
    primaryVars = {'pressure', 'sW', 'sulfate', 'calcium','strontium','barium'};
end

% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations), polymer concentration and well rates +
% bhp.

% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;
sat  = {sW, sO};
sat0 = {sW0, sO0};


% Update state with AD-variables
state = model.setProps(state  , {'s', 'pressure', 'sulfate', 'calcium', 'strontium', 'barium', 'sodium', 'magnesium', 'carbonate','chlorine'},...
                                {sat , p , c_SO4, c_Ca, c_Sr, c_Ba, c_Na, c_Mg, c_CO3, c_Cl});
state0 = model.setProps(state0, {'s', 'pressure', 'sulfate', 'calcium', 'strontium', 'barium', 'sodium', 'magnesium', 'carbonate','chlorine'},...
                                {sat0, p0, c0_SO4, c0_Ca, c0_Sr, c0_Ba, c0_Na, c0_Mg, c0_CO3, c0_Cl});
% Set up properties
% state = model.initStateFunctionContainers(state);
state = model.initStateFunctionContainers(state);

% Compute transmissibility
T = op.T;

% Gravity contribution
gdz = model.getGravityGradient();

[b, pv] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume');
[b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');
[phaseFlux, flags] = model.getProps(state, 'PhaseFlux',  'PhaseUpwindFlag');
[pressures, mob, rho] = model.getProps(state, 'PhasePressures', 'Mobility', 'Density');
bO=deal(b{1,1});
bW=deal(b{1,2});
bO0=deal(b0{1,1});
bW0=deal(b0{1,2});

[vW, vO] = deal(phaseFlux{:});
[upcw, upco] = deal(flags{:});

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, []);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end


% Set up properties
% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.

bOvO = op.faceUpstr(upco, bO).*vO;

%Water fraction
% cw  = 1-c_Ca-c_Ba-c_Sr-c_Na-c_Mg-c_CO3-c_SO4-c_Cl;
% cw0 = 1-c0_Ca-c0_Ba-c0_Sr-c0_Na-c0_Mg-c0_CO3-c0_SO4-c0_Cl;
% bWvW = op.faceUpstr(upcw, cw).*vW;
bWvW = op.faceUpstr(upcw, bW).*vW;

vSO4 = op.faceUpstr(upcw, c_SO4).*vW;
bWvSO4 = op.faceUpstr(upcw, bW).*vSO4;

vCa = op.faceUpstr(upcw, c_Ca).*vW;
bWvCa = op.faceUpstr(upcw, bW).*vCa;

vSr = op.faceUpstr(upcw, c_Sr).*vW;
bWvSr = op.faceUpstr(upcw, bW).*vSr;

vBa = op.faceUpstr(upcw, c_Ba).*vW;
bWvBa = op.faceUpstr(upcw, bW).*vBa;

vNa = op.faceUpstr(upcw, c_Na).*vW;
bWvNa = op.faceUpstr(upcw, bW).*vNa;

vMg = op.faceUpstr(upcw, c_Mg).*vW;
bWvMg = op.faceUpstr(upcw, bW).*vMg;

vCO3 = op.faceUpstr(upcw, c_CO3).*vW;
bWvCO3 = op.faceUpstr(upcw, bW).*vCO3;

vCl = op.faceUpstr(upcw, c_Cl).*vW;
bWvCl = op.faceUpstr(upcw, bW).*vCl;

% Conservation of mass for water
alpha = model.pvCorrector;
water = (op.pv/dt).*( bW.*sW - alpha.*bW0.*sW0 ) + op.Div(bWvW);

% Conservation of mass for oil
oil = (op.pv/dt).*( bO.*sO - alpha.*bO0.*sO0 ) + op.Div(bOvO);

% Conservation of ions in water
D_SO4 = 1e-11; % Diffusion coefficient (m^2/s)
D_Ca  = 1e-11; % Diffusion coefficient (m^2/s)
D_Sr  = 1e-11; % Diffusion coefficient (m^2/s)
D_Ba  = 1e-11; % Diffusion coefficient (m^2/s)
D_Na  = 1e-11; % Diffusion coefficient (m^2/s)
D_Mg  = 1e-11; % Diffusion coefficient (m^2/s)
D_CO3 = 1e-11; % Diffusion coefficient (m^2/s)
D_Cl  = 1e-11; % Diffusion coefficient (m^2/s)

diff_SO4 = op.Div(D_SO4.*op.faceUpstr(upcw, bW).*op.Grad(c_SO4)); % diffusion term
diff_Ca  = op.Div(D_Ca.*op.faceUpstr(upcw, bW).*op.Grad(c_Ca));   % diffusion term 
diff_Sr  = op.Div(D_Sr.*op.faceUpstr(upcw, bW).*op.Grad(c_Sr));   % diffusion term
diff_Ba  = op.Div(D_Ba.*op.faceUpstr(upcw, bW).*op.Grad(c_Ba));   % diffusion term
diff_Na  = op.Div(D_Na.*op.faceUpstr(upcw, bW).*op.Grad(c_Na));   % diffusion term
diff_Mg  = op.Div(D_Mg.*op.faceUpstr(upcw, bW).*op.Grad(c_Mg));   % diffusion term
diff_CO3 = op.Div(D_CO3.*op.faceUpstr(upcw, bW).*op.Grad(c_CO3));   % diffusion term
diff_Cl  = op.Div(D_Cl.*op.faceUpstr(upcw, bW).*op.Grad(c_Cl));   % diffusion term


sulfate = (op.pv/dt).*(bW.*sW.*c_SO4-alpha.*bW0.*sW0.*c0_SO4) + op.Div(bWvSO4) - diff_SO4;
calcium = (op.pv/dt).*(bW.*sW.*c_Ca-alpha.*bW0.*sW0.*c0_Ca) + op.Div(bWvCa) - diff_Ca ;
strontium = (op.pv/dt).*(bW.*sW.*c_Sr-alpha.*bW0.*sW0.*c0_Sr) + op.Div(bWvSr) - diff_Sr;
barium = (op.pv/dt).*(bW.*sW.*c_Ba-alpha.*bW0.*sW0.*c0_Ba) + op.Div(bWvBa) - diff_Ba;
sodium = (op.pv/dt).*(bW.*sW.*c_Na-alpha.*bW0.*sW0.*c0_Na) + op.Div(bWvNa) - diff_Na;
magnesium = (op.pv/dt).*(bW.*sW.*c_Mg-alpha.*bW0.*sW0.*c0_Mg) + op.Div(bWvMg) - diff_Mg;
carbonate = (op.pv/dt).*(bW.*sW.*c_CO3-alpha.*bW0.*sW0.*c0_CO3) + op.Div(bWvCO3) - diff_CO3;
chlorine = (op.pv/dt).*(bW.*sW.*c_Cl-alpha.*bW0.*sW0.*c0_Cl) + op.Div(bWvCl) - diff_Cl;



if ~opt.resOnly
    epsilon = 1.e-8;
    % the first way is based on the diagonal values of the resulting
    % Jacobian matrix
    dSO4 = diag(sulfate.jac{3});
    dCa = diag(calcium.jac{4});
    dSr = diag(strontium.jac{5});
    dBa = diag(barium.jac{6});
    dNa = diag(sodium.jac{7});
    dMg = diag(magnesium.jac{8});
    dCO3 = diag(carbonate.jac{9});
    dCl = diag(chlorine.jac{10});
    epsSO4 = sqrt(epsilon)*mean(abs(dSO4));
    epsCa = sqrt(epsilon)*mean(abs(dCa));
    epsSr = sqrt(epsilon)*mean(abs(dSr));
    epsBa = sqrt(epsilon)*mean(abs(dBa));
    epsNa = sqrt(epsilon)*mean(abs(dNa));
    epsMg = sqrt(epsilon)*mean(abs(dMg));
    epsCO3 = sqrt(epsilon)*mean(abs(dCO3));
    epsCl = sqrt(epsilon)*mean(abs(dCl));
    % sometimes there is no water in the whole domain
    if (epsSO4 == 0.)
        epsSO4 = epsilon;
    end
    if (epsCa == 0.)
        epsCa = epsilon;
    end
    if (epsSr == 0.)
        epsSr = epsilon;
    end
    if (epsBa == 0.)
        epsBa = epsilon;
    end
    if (epsNa == 0.)
        epsNa = epsilon;
    end
    if (epsMg == 0.)
        epsMg = epsilon;
    end
    if (epsCO3 == 0.)
        epsCO3 = epsilon;
    end
    if (epsCl == 0.)
        epsCl = epsilon;
    end
    % bad marks the cells prolematic in evaluating Jacobian
    badSO4 = abs(dSO4) < epsSO4;
    badCa = abs(dCa) < epsCa;
    badSr = abs(dSr) < epsSr;
    badBa = abs(dBa) < epsBa;
    badNa = abs(dNa) < epsNa;
    badMg = abs(dMg) < epsMg;
    badCO3 = abs(dCO3) < epsCO3;
    badCl = abs(dCl) < epsCl;
    % the other way is to choose based on the water saturation
    sulfate(badSO4) = c_SO4(badSO4);
    calcium(badCa) = c_Ca(badCa);
    strontium(badSr) = c_Sr(badSr);
    barium(badBa) = c_Ba(badBa);
    sodium(badNa) = c_Na(badNa);
    magnesium(badMg) = c_Mg(badMg);
    carbonate(badCO3) = c_CO3(badCO3);
    chlorine(badCl) = c_Cl(badCl);
end
eqs   = {water, oil ,sulfate, calcium, strontium, barium, sodium, magnesium, carbonate, chlorine};
names = {'water', 'oil','sulfate', 'calcium', 'strontium', 'barium', 'sodium', 'magnesium', 'carbonate', 'chlorine'};
types = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
components = {c_SO4, c_Ca, c_Sr, c_Ba, c_Na, c_Mg, c_CO3, c_Cl};
dissolved = {};

sat = {sW, sO};
[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 pressures, sat, mob, rho, ...
                                                                 dissolved, components, ...
                                                                 drivingForces);

% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, components, dt, opt);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end