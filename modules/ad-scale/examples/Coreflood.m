mrstModule add ad-core ad-blackoil ad-scale ad-props

% Define grid
[NX,NY,NZ]=deal(10,1,1);
[dx,dy,dz]=deal(0.153,0.0337,0.0337);
G=computeGeometry(cartGrid([NX,NY,NZ],[dx,dy,dz]));
initial_k = 9*milli*darcy;
initial_phi = 0.10;

%Making rock and fluid
rock = makeRock(G, initial_k, initial_phi);

fluid = initSimpleADIFluid('mu',    [0.5, 4]*centi*poise, ...
                           'rho',   [1100, 800]*kilogram/meter^3, ...
                           'cR',    1e-6/psia, ...
                           'c',     [1e-6/psia 4e-6/psia],...
                           'phases', 'wo');
                       

% %Rel_Perm Table
Sw_kr = [0 0.2688 0.31027658 0.39909346 0.46871307 0.52582234 0.56992975 0.6080874 0.64107817 0.66770833 0.68900712 0.70263715 0.7145 0.7312 1];
krw_table = [0 0 0.03953 0.094104 0.12838 0.153072 0.171206 0.187268 0.202616 0.2174 0.232489 0.245156 0.26 0.265 0.265];
kro_table = [0.45 0.45 0.301756 0.18166 0.113294 0.070646 0.044936 0.027383 0.015578 0.008318 0.004015 0.002006 0.000774 0.000001 0];
 
fluid.krW = @(sw) interpTable(Sw_kr,krw_table,sw);
fluid.krO = @(so) interpTable(Sw_kr,kro_table,1-so); 
% %Pc implement
SW_table = [0 0.2688 0.291085 0.31337 0.335655 0.3579 0.380225 0.40251 0.424795 0.44708 0.469365 0.49165 0.513935 0.53622 0.558505 0.58079 0.603075 0.62536 0.647645 0.66993 0.692215 0.701129 0.710043 0.7145];
pc_table = psia.*[0 0 -0.3067 -0.42101 -0.57792 -0.79331 -1.08897 -1.49484 -2.05196 -2.81674 -3.86654 -5.30761 -7.28577 -10.0012 -13.7287 -18.8454 -25.8691 -35.5105 -48.7454 -66.9129 -91.8515 -104.26 -118.344 -136.4];
fluid.pcOW = @(sw)interpTable(SW_table,pc_table,sw); 
                      
% Construct reservoir model
gravity reset on
sWi = 0.2688;
sOi = 1-sWi;
initial_Pressure = 100*psia;
state0         = initResSol(G, initial_Pressure, [sWi, 1-sWi]);
ppm0_Ca  = 9500;
ppm0_Sr  = 200;
ppm0_Ba  = 5;
ppm0_SO4 = 100;
ppm0_Na  = 70015;
ppm0_Mg  = 20;
ppm0_CO3 = 10;
ppm0_Cl  = 129022;
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
model.toleranceMB=1e-4;
model.toleranceCNV=5e-3;

%schedule
simTime = 3600;
nstep   =1000 ;
% startsteps_dt = 2;
% startsteps = repmat(startsteps_dt, 200, 1);
% reststeps_dt = 2;
% reststeps = repmat(reststeps_dt, nstep-200, 1);
% timesteps = [startsteps;reststeps];
dt = 200;
timesteps = repmat(dt, nstep, 1);

%Boundary conditions
bc1 = fluxside([], G, 'Left', 3.12768*centi^3/hour, 'sat', [1 0]);
bc2 = fluxside([], G, 'Left', 3.12768*centi^3/hour, 'sat', [1 0]);
ppm_Ca_injection1  = 9500;
ppm_Sr_injection1  = 200;
ppm_Ba_injection1  = 5;
ppm_SO4_injection1 = 100;
ppm_Na_injection1  = 70015;
ppm_Mg_injection1  = 20;
ppm_CO3_injection1 = 10;
ppm_Cl_injection1  = 129022;
sigma_injection1 = 1e-6*(ppm_Ca_injection1+ppm_Sr_injection1+ppm_Ba_injection1+ppm_SO4_injection1+ppm_Na_injection1+ppm_Mg_injection1+ppm_CO3_injection1+ppm_Cl_injection1);
bc1.c_SO4=(ppm_SO4_injection1*1e-6/(1-sigma_injection1)).*ones(size(bc1.sat,1), 1);
bc1.c_Ca =(ppm_Ca_injection1*1e-6/(1-sigma_injection1)).*ones(size(bc1.sat,1), 1);
bc1.c_Sr =(ppm_Sr_injection1*1e-6/(1-sigma_injection1)).*ones(size(bc1.sat,1), 1);
bc1.c_Ba =(ppm_Ba_injection1*1e-6/(1-sigma_injection1))*ones(size(bc1.sat,1), 1);
bc1.c_Na =(ppm_Na_injection1*1e-6/(1-sigma_injection1)).*ones(size(bc1.sat,1), 1);
bc1.c_Mg =(ppm_Mg_injection1*1e-6/(1-sigma_injection1)).*ones(size(bc1.sat,1), 1);
bc1.c_CO3=(ppm_CO3_injection1*1e-6/(1-sigma_injection1)).*ones(size(bc1.sat,1), 1);
bc1.c_Cl =(ppm_Cl_injection1*1e-6/(1-sigma_injection1)).*ones(size(bc1.sat,1), 1);

ppm_Ca_injection2  = 75;
ppm_Sr_injection2  = 0;
ppm_Ba_injection2  = 0;
ppm_SO4_injection2 = 13240;
ppm_Na_injection2  = 9720;
ppm_Mg_injection2  = 1306;
ppm_CO3_injection2 = 200;
ppm_Cl_injection2  = 20955;
sigma_injection2 = 1e-6*(ppm_Ca_injection2+ppm_Sr_injection2+ppm_Ba_injection2+ppm_SO4_injection2+ppm_Na_injection2+ppm_Mg_injection2+ppm_CO3_injection2+ppm_Cl_injection2);
bc2.c_SO4=(ppm_SO4_injection2*1e-6/(1-sigma_injection2)).*ones(size(bc2.sat,1), 1);
bc2.c_Ca =(ppm_Ca_injection2*1e-6/(1-sigma_injection2)).*ones(size(bc2.sat,1), 1);
bc2.c_Sr =(ppm_Sr_injection2*1e-6/(1-sigma_injection2)).*ones(size(bc2.sat,1), 1);
bc2.c_Ba =(ppm_Ba_injection2*1e-6/(1-sigma_injection2))*ones(size(bc2.sat,1), 1);
bc2.c_Na =(ppm_Na_injection2*1e-6/(1-sigma_injection2)).*ones(size(bc2.sat,1), 1);
bc2.c_Mg =(ppm_Mg_injection2*1e-6/(1-sigma_injection2)).*ones(size(bc2.sat,1), 1);
bc2.c_CO3=(ppm_CO3_injection2*1e-6/(1-sigma_injection2)).*ones(size(bc2.sat,1), 1);
bc2.c_Cl =(ppm_Cl_injection2*1e-6/(1-sigma_injection2)).*ones(size(bc2.sat,1), 1);

bc1 = pside(bc1, G, 'Right', 100*psia, 'sat', [0 1]);
bc2 = pside(bc2, G, 'Right', 100*psia, 'sat', [0 1]);

schedule = simpleSchedule(timesteps);
schedule.step.control(201:1000)=2;
tmp = cell(2,1);
schedule.control=struct('W',tmp,'bc',tmp,'src',tmp);
schedule.control(1).bc=bc1;
schedule.control(2).bc=bc2;


fn = getPlotAfterStep(state0, model, schedule, ...
                      'plotReservoir', true, 'view', [20, 8], ...
                      'field', 'c_SO4');

[wellSols, states, rockProps,report] = simulateScheduleAD(state0, model, schedule,'afterStepFn', fn);

CoreData = CoreDataProcess(states, nstep,NX,dx,rockProps);
ions(1:nstep,1) = CoreData.Na_out;
ions(1:nstep,2) = CoreData.Cl_out;
ions(1:nstep,3) = CoreData.SO4_out;
ions(1:nstep,4) = CoreData.Ca_out;
ions(1:nstep,5) = CoreData.Mg_out;
ions(1:nstep,6) = CoreData.CO3_out;
ions(1:nstep,7) = CoreData.Sr_out;
ions(1:nstep,8) = CoreData.Ba_out;

T = [0:200:nstep*dt]/3600;

figure
semilogy(T,ions(:,1))
grid
hold
semilogy(T,ions(:,2))
semilogy(T,ions(:,3))
semilogy(T,ions(:,4))
semilogy(T,ions(:,5))
semilogy(T,ions(:,6))
semilogy(T,ions(:,7))
semilogy(T,ions(:,8))
legend('Na','Cl','SO4','Ca','Mg','CO3','Sr','Ba')
xlabel('Time (hour)')
ylabel('Outlet Concentration')

