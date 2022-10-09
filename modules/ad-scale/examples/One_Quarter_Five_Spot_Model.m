%% One_Quarte_Five_Spot
% 3000 days formation brine injection followed by 3000 days seawater
% injectio
mrstModule add ad-core ad-blackoil ad-scale ad-props
% Define grid
[NX,NY,NZ]=deal(10,10,1);
[dx,dy,dz]=deal(141.42,141.42,15);
G=computeGeometry(cartGrid([NX,NY,NZ],[dx,dy,dz]));
initial_k = [9,9,0.9].*milli*darcy;
initial_phi = 0.1;

%Making rock and fluid
rock = makeRock(G, initial_k, initial_phi);

fluid = initSimpleADIFluid('mu',    [1, 4]*centi*poise, ...
                           'rho',   [1200, 800]*kilogram/meter^3, ...
                           'cR',    1e-6/psia, ...
                           'c',     [1e-6/psia 6e-6/psia],...
                           'phases', 'wo');

%Rel_Perm Table
Sw_kr = [-1 0.2688 0.31027658 0.39909346 0.46871307 0.52582234 0.56992975 0.6080874 0.64107817 0.66770833 0.68900712 0.70263715 0.7145 0.7312 1];
krw_table = [0 0 0.03953 0.094104 0.12838 0.153072 0.171206 0.187268 0.202616 0.2174 0.232489 0.245156 0.26 0.265 0.265];
kro_table = [0.45 0.45 0.301756 0.18166 0.113294 0.070646 0.044936 0.027383 0.015578 0.008318 0.004015 0.002006 0.000774 0 0];
 
fluid.krW = @(sw) interpTable(Sw_kr,krw_table,sw);
fluid.krO = @(so) interpTable(Sw_kr,kro_table,1-so); 
%Pc implement
SW_table = [-1 0.2688 0.291085 0.31337 0.335655 0.3579 0.380225 0.40251 0.424795 0.44708 0.469365 0.49165 0.513935 0.53622 0.558505 0.58079 0.603075 0.62536 0.647645 0.66993 0.692215 0.701129 0.710043 0.7145];
pc_table = psia.*[0 0 -0.3067 -0.42101 -0.57792 -0.79331 -1.08897 -1.49484 -2.05196 -2.81674 -3.86654 -5.30761 -7.28577 -10.0012 -13.7287 -18.8454 -25.8691 -35.5105 -48.7454 -66.9129 -91.8515 -104.26 -118.344 -136.4];
fluid.pcOW = @(sw)interpTable(SW_table,pc_table,sw); 

% Construct reservoir model
gravity reset on
sWi = 0.3;
initial_Pressure = 4980*psia;
state0         = initResSol(G, initial_Pressure, [sWi, 1-sWi]);
ppm_Ca  = 9500;
ppm_Sr  = 200;
ppm_Ba  = 5;
ppm_SO4 = 100;
ppm_Na  = 70015;
ppm_Mg  = 10;
ppm_CO3 = 5;
ppm_Cl  = 129022;
sigma_ions = (ppm_Ca + ppm_Sr + ppm_Ba + ppm_SO4 + ppm_Na + ppm_Mg + ppm_CO3 + ppm_Cl)*1e-6;
state0.c_SO4   = (ppm_SO4*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Ca    = (ppm_Ca*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Sr    = (ppm_Sr*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Ba    = (ppm_Ba*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Na    = (ppm_Na*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Mg    = (ppm_Mg*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_CO3   = (ppm_CO3*1e-6/(1-sigma_ions))*ones(G.cells.num,1);
state0.c_Cl    = (ppm_Cl*1e-6/(1-sigma_ions))*ones(G.cells.num,1);

model = OilWaterScaleModel(G, rock, fluid);
model.initial_props.poro = initial_phi;
model.initial_props.sWi = sWi;
model.initial_props.T = 120; %temperature
model.initial_props.P = initial_Pressure;
model.initial_props.perm = initial_k;
model.pvCorrector = ones(NX*NY*NZ,1);

% model.gravity=-model.gravity;
if NZ>1
    nstep0 = 1;
    dt = 10000*day;
    timesteps = repmat(dt, nstep0, 1);
    schedule0 = simpleSchedule(timesteps);
    [~, state_initial,~, ~] = simulateScheduleAD(state0, model, schedule0);
    state0.pressure  = state_initial{nstep0,1}.pressure;
    state0.flux      = state_initial{nstep0,1}.flux;
    state0.s         = state_initial{nstep0,1}.s;
    state0.c_Ca      = state_initial{nstep0,1}.c_Ca;
    state0.c_Mg      = state_initial{nstep0,1}.c_Mg;
    state0.c_Sr      = state_initial{nstep0,1}.c_Sr;
    state0.c_Ba      = state_initial{nstep0,1}.c_Ba;
    state0.c_SO4     = state_initial{nstep0,1}.c_SO4;
    state0.c_Cl      = state_initial{nstep0,1}.c_Cl;
    state0.c_Na      = state_initial{nstep0,1}.c_Na;
    state0.c_CO3     = state_initial{nstep0,1}.c_CO3;
    clear dt nstep0 timesteps schedule0
    clc
end

simTime = 6000*day;
nstep   = 200;
startsteps_dt = 30*day;
startsteps = repmat(startsteps_dt, 100, 1);
reststeps_dt = 30*day;
reststeps = repmat(reststeps_dt, nstep-100, 1);
timesteps = [startsteps;reststeps];


injRate = 100*stb/day;


W1 = verticalWell([],G,rock,NX,NY,1:NZ,...
                  'Type', 'bhp', 'Val', 3000*psia,...
                  'Name','Producer','comp_i',[0.5 0.5], 'sign', -1);

W1 = verticalWell(W1,G,rock,1,1,1:NZ,...
                  'Type', 'bhp', 'Val', 5000*psia,...
                  'Name','Injector','comp_i',[1 0], 'sign', 1);


W1(1).c_SO4  = 0;
W1(1).c_Ca   = 0;
W1(1).c_Sr   = 0;
W1(1).c_Ba   = 0;
W1(1).c_Na   = 0;
W1(1).c_Mg   = 0;
W1(1).c_CO3  = 0;
W1(1).c_Cl   = 0;
W1(2).c_SO4  = (ppm_SO4*1e-6/(1-sigma_ions));
W1(2).c_Ca   = (ppm_Ca*1e-6/(1-sigma_ions));
W1(2).c_Sr   = (ppm_Sr*1e-6/(1-sigma_ions));
W1(2).c_Ba   = (ppm_Ba*1e-6/(1-sigma_ions));
W1(2).c_Na   = (ppm_Na*1e-6/(1-sigma_ions));
W1(2).c_Mg   = (ppm_Mg*1e-6/(1-sigma_ions));
W1(2).c_CO3  = (ppm_CO3*1e-6/(1-sigma_ions));
W1(2).c_Cl   = (ppm_Cl*1e-6/(1-sigma_ions));

W2 = verticalWell([],G,rock,NX,NY,1:NZ,...
                  'Type', 'bhp', 'Val', 3000*psia,...
                  'Name','Producer','comp_i',[0.5 0.5], 'sign', -1);

W2 = verticalWell(W2,G,rock,1,1,1:NZ,...
                  'Type', 'bhp', 'Val', 5000*psia,...
                  'Name','Injector','comp_i',[1 0], 'sign', 1);
ppm_SO4_injection = 13240;
ppm_Ca_injection  = 75;
ppm_Sr_injection  = 0;
ppm_Ba_injection  = 0;
ppm_Na_injection  = 9720;
ppm_Mg_injection  = 1306;
ppm_CO3_injection = 200;
ppm_Cl_injection  = 20955;
sigma_injection = 1e-6*(ppm_Ca_injection+ppm_Sr_injection+ppm_Ba_injection+ppm_SO4_injection+ppm_Na_injection+ppm_Mg_injection+ppm_CO3_injection+ppm_Cl_injection);
W2(1).c_SO4  = 0;
W2(1).c_Ca   = 0;
W2(1).c_Sr   = 0;
W2(1).c_Ba   = 0;
W2(1).c_Na   = 0;
W2(1).c_Mg   = 0;
W2(1).c_CO3  = 0;
W2(1).c_Cl   = 0;
W2(2).c_SO4  = (ppm_SO4_injection*1e-6/(1-sigma_injection));
W2(2).c_Ca   = (ppm_Ca_injection*1e-6/(1-sigma_injection));
W2(2).c_Sr   = (ppm_Sr_injection*1e-6/(1-sigma_injection));
W2(2).c_Ba   = (ppm_Ba_injection*1e-6/(1-sigma_injection));
W2(2).c_Na   = (ppm_Na_injection*1e-6/(1-sigma_injection));
W2(2).c_Mg   = (ppm_Mg_injection*1e-6/(1-sigma_injection));
W2(2).c_CO3  = (ppm_CO3_injection*1e-6/(1-sigma_injection));
W2(2).c_Cl   = (ppm_Cl_injection*1e-6/(1-sigma_injection));
          
              
% state0 =Initializer(model,state0);

schedule = simpleSchedule(timesteps);
tmp = cell(2,1);
schedule.step.control(100:600)=2;
schedule.control=struct('W',tmp,'bc',tmp,'src',tmp);
schedule.control(1).W=W1;
schedule.control(2).W=W2;
% [wellSols, states,rockProps, report] = simulateScheduleAD(state0, model, schedule);
fn = getPlotAfterStep(state0, model, schedule, ...
                      'plotReservoir', true, ...
                      'field', 'c_SO4',...
                      'wells',W1,'view',[20 50]);

[wellSols, states, rockProps] = simulateScheduleAD(state0, model, schedule,'afterStepFn', fn);
                                    
                                    
                                    
                                    
wellData = WellDataProcess(wellSols,nstep,reststeps_dt);
ReservoirData = ReservoirDataProcess(states,nstep);
% qOt   = deal(zeros(nstep,1));
% p_ave = deal(zeros(nstep,1));
% phi   = deal(zeros(nstep,1));
% perm  = deal(cell(nstep,1));
% PV = sum(G.cells.volumes)*initial_phi;
% PV = convertTo(PV,stb);                                    
 
                                    
t = 1:nstep;
% plot(t,convertTo(qOr,(stb/day)),'*')
% plot(t,convertTo(qWr,(stb/day)),'-')
% xlabel('timestep')
% ylabel('Rate (stb/day)')
% legend('Oil Production Rate','Water Production Rate')  

figure
grid
hold
plot(t,wellData.Ca,'b*')
plot(t,wellData.Sr,'r-')
plot(t,wellData.Ba,'y+')
plot(t,wellData.SO4,'k.-')
plot(t,wellData.Na,'gx')
plot(t,wellData.Cl,'ms')
plot(t,wellData.CO3,'b--')
plot(t,wellData.Mg,'cd')
xlabel('time step')
ylabel('ppm')
legend('Ca','Sr','Ba','SO4','Na','Cl','CO3','Mg')                               