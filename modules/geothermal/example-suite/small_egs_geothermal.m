function [description, options, state0, model, schedule, plotOptions] = small_egs_geothermal(varargin)
    description = 'Small, enhanced geothermal system';
    options = struct('ncells', 6, 'nstep', 30, 'nlayers', 8);
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Module dependencies
    require ad-core ad-props ad-blackoil geothermal compositional upr
    gravity reset on
    rng(20200629)
    % Load fracture dataset
    pth  = fullfile(getDatasetPath('geothermal'), 'egs-fractures.mat');
    data = load(pth);
    xmax = [50,50,15];                    % Physical dimensions
    dx   = min(xmax(1:2))/options.ncells; % Approximate cell size
    G2D = pebiGrid2D(dx, xmax(1:2), ...
            'cellConstraints', data.fractures, ... % Fractures
            'CCRefinement'   , true          , ... % Refine fractures
            'CCFactor'       , 0.1           );    % Rel fracture cell size
    % Make layered grid
    layers = diff(linspace(0, xmax(3), options.nlayers + 1));
    G = makeLayeredGrid(G2D, layers);
    G = computeGeometry(G);
    % Repeat fracture tags
    G.cells.tag = repmat(G.cells.tag, options.nlayers, 1);
    % Place wells in opposite ends of fracture newtork
    fracture_cells = G.cells.tag;
    xf = sum(G.cells.centroids(fracture_cells,1:2).^2,2);
    [~, xf] = sort(xf);
    inj  = xf(1:options.nlayers);
    prod = xf(end-options.nlayers+1:end);
    map = find(fracture_cells);
    inj = map(inj); prod = map(prod);
    % Make rock
    rock = makeRock(G, 1*milli*darcy, 0.05);
    rock.perm(fracture_cells) = 1e4*milli*darcy; % Fracture permeability
    rock.poro(fracture_cells) = 0.8;             % Fracture porosity
    % Add thermal properties
    Watt = joule/second;
    rock = addThermalRockProps(rock, 'lambdaR', 2*Watt/(meter*Kelvin)       , ...
                                     'rhoR'   , 2700*kilogram/meter^3       , ...
                                     'CpR'    , 1000*joule/(kilogram*Kelvin));
    K0   = 273.15*Kelvin;
    % Make fluid
    fluid = initSimpleADIFluid(                     ...
                   'phases', 'W'                  , ... % Water only
                   'mu'    , 0.5*centi*poise      , ... % Viscosity
                   'rho'   , 1000*kilogram/meter^3, ... % Reference density
                   'c'     , 4.4e-10/Pascal       , ... % Compressibility
                   'pRef'  , 1*atm                );    % Ref pressure
    % Add thermal properties
    fluid = addThermalFluidProps(fluid            , ...
                'Cp'     , 4.2*joule/(gram*Kelvin), ... % Heat capacity
                'lambdaF', 0.6*Watt/(meter*Kelvin), ... % Thermal cond
                'cT'     , 207e-6/Kelvin          , ... % Thermal expansion
                'TRef'   , K0 + 10*Kelvin         );    % Reference temp
    % Make model
    model = GeothermalModel(G, rock, fluid);
    model.extraStateOutput = true;
    model.outputFluxes     = true;
    % Initial state
    K0   = 273.15*Kelvin;
    Tres = K0 + 95*Kelvin;
    Tinj = K0 + 10*Kelvin;
    state0   = initResSol(G, 300*barsa, 1);
    state0.T = repmat(Tres, G.cells.num, 1);
    % Inject 1 fracture pore volume per year for 50 years
    time = 25*year;
    pv   = poreVolume(G, rock);
    rate = 0.5*sum(pv(fracture_cells))/year;
    % Set up wells
    W = addWell([], G, rock, inj, ...
                'type', 'rate', 'val', rate, 'compi', 1, 'name', 'Inj');
    W = addWell(W, G, rock, prod, ...
                'type', 'bhp', 'val', 50*barsa, 'compi', 1, 'name', 'Prod');
    W = addThermalWellProps(W, 'T', Tinj);
    % Make schedule
    schedule = simpleSchedule(rampupTimesteps(time, time/options.nstep), 'W', W);
    % Plotting
    plotOptions = {'View', [65, 30], ...
                   'Size', [700, 500], ...
                   'PlotBoxAspectRatio', [1, xmax(2)/xmax(1), 15/xmax(1)]};
end