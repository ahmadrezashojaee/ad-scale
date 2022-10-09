mrstModule add ad-core ad-blackoil ad-scale ad-props

% Define grid
[NX,NY,NZ]=deal(10,1,1);
[dx,dy,dz]=deal(0.153,0.0337,0.0337);
G=computeGeometry(cartGrid([NX,NY,NZ],[dx,dy,dz]));
initial_k = 25*milli*darcy;
initial_phi = 0.17;

%Making rock and fluid
rock = makeRock(G, initial_k, initial_phi);

fluid = initSimpleADIFluid('mu',    [0.3, 4]*centi*poise, ...
                           'rho',   [1100, 800]*kilogram/meter^3, ...
                           'cR',    1e-6/psia, ...
                           'c',     [1e-6/psia 12e-6/psia],...
                           'phases', 'wo',...
                           'smin',[0 0],...
                           'n',[0 2]);
                      
% Construct reservoir model
gravity reset on
sWi = 1;
sOi = 1-sWi;
initial_Pressure = 55*psia;
state0         = initResSol(G, initial_Pressure, [sWi, 1-sWi]);
ppm0_Ca  = 11236;
ppm0_Sr  = 0;
ppm0_Ba  = 0;
ppm0_SO4 = 215;
ppm0_Na  = 40288;
ppm0_Mg  = 2815;
ppm0_CO3 = 0;
ppm0_Cl  = 107309;
sigma_ions = (ppm0_Ca + ppm0_Sr + ppm0_Ba + ppm0_SO4 + ppm0_Na + ppm0_Mg + ppm0_CO3 + ppm0_Cl)*1e-6;
state0.c_SO4   = (ppm0_SO4*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Ca    = (ppm0_Ca*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Sr    = (ppm0_Sr*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Ba    = (ppm0_Ba*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Na    = (ppm0_Na*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Mg    = (ppm0_Mg*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_CO3   = (ppm0_CO3*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Cl    = (ppm0_Cl*1e-6/(1-sigma_ions))*ones(G.cells.num,1);

model = OilWaterScaleModel(G, rock, fluid);
model.initial_props.poro = initial_phi;
model.initial_props.sWi = sWi;
model.initial_props.T = 120; %temperature
model.initial_props.P = initial_Pressure;
model.initial_props.perm = initial_k;
model.pvCorrector = ones(NX*NY*NZ,1);
model.toleranceMB=1e-5;
model.toleranceCNV=1e-3;

%schedule
simTime = 3600;
nstep   =300 ;
% startsteps_dt = 2;
% startsteps = repmat(startsteps_dt, 200, 1);
% reststeps_dt = 2;
% reststeps = repmat(reststeps_dt, nstep-200, 1);
% timesteps = [startsteps;reststeps];
dt = 118.5;
timesteps = repmat(dt, nstep, 1);

%Boundary conditions

bc = fluxside([], G, 'Left', 6*centi^3/hour, 'sat', [1 0]);
ppm_Ca_injection  = 507;
ppm_Sr_injection  = 0;
ppm_Ba_injection  = 0;
ppm_SO4_injection = 3485;
ppm_Na_injection  = 12892;
ppm_Mg_injection  = 1519;
ppm_CO3_injection = 0;
ppm_Cl_injection  = 26578;
sigma_injection = 1e-6*(ppm_Ca_injection+ppm_Sr_injection+ppm_Ba_injection+ppm_SO4_injection+ppm_Na_injection+ppm_Mg_injection+ppm_CO3_injection+ppm_Cl_injection);
bc.c_SO4=(ppm_SO4_injection*1e-6/(1-sigma_injection)).*ones(size(bc.sat,1), 1);
bc.c_Ca =(ppm_Ca_injection*1e-6/(1-sigma_injection)).*ones(size(bc.sat,1), 1);
bc.c_Sr =(ppm_Sr_injection*1e-6/(1-sigma_injection)).*ones(size(bc.sat,1), 1);
bc.c_Ba =(ppm_Ba_injection*1e-6/(1-sigma_injection))*ones(size(bc.sat,1), 1);
bc.c_Na =(ppm_Na_injection*1e-6/(1-sigma_injection)).*ones(size(bc.sat,1), 1);
bc.c_Mg =(ppm_Mg_injection*1e-6/(1-sigma_injection)).*ones(size(bc.sat,1), 1);
bc.c_CO3=(ppm_CO3_injection*1e-6/(1-sigma_injection)).*ones(size(bc.sat,1), 1);
bc.c_Cl =(ppm_Cl_injection*1e-6/(1-sigma_injection)).*ones(size(bc.sat,1), 1);

bc = pside(bc, G, 'Right', 55*psia, 'sat', [0 1]);


schedule = simpleSchedule(timesteps, 'bc', bc);
fn = getPlotAfterStep(state0, model, schedule, ...
                      'plotReservoir', true, 'view', [20, 8], ...
                      'field', 'c_SO4');

[wellSols, states, rockProps,Reports] = simulateScheduleAD(state0, model, schedule,'afterStepFn', fn);

pv = [1:300]/150;
[CoreData] = CoreDataProcess(states, nstep,NX,dx,rockProps);
benchmark1;
plotting;