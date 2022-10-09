%Five Spot Model
mrstModule add ad-core ad-blackoil ad-scale ad-props

% Define grid
[NX,NY,NZ]=deal(11,11,1);
[dx,dy,dz]=deal(600,600,5);
G=computeGeometry(cartGrid([NX,NY,NZ],[dx,dy,dz]));
initial_k = 100*milli*darcy;
initial_phi = 0.20;

rock = makeRock(G, initial_k, initial_phi);


fluid = initSimpleADIFluid('mu',    [1, 4]*centi*poise, ...
                           'rho',   [1100, 800]*kilogram/meter^3, ...
                           'n',     [3, 2], ...
                           'cR',    1e-7/psia, ...
                           'c',     [1e-6*psia^-1 12e-6*psia^-1],...
                           'smin',  [0.2 0.2],...
                           'phases', 'wo');
%Pc implement
SW_table = [-0.1 0:0.05:1 2];
pc_table = [500 500 300 200 150 10 5 3.5 3 2 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 0.2 -0.1 -0.2 -0.3 -0.3];
fluid.pcOW = @(sw)interpTable(SW_table,pc_table,sw);                     
% Construct reservoir model
gravity reset on
sWi = 0.2;
initial_Pressure = 3000*psia;
state0         = initResSol(G, initial_Pressure, [sWi, 1-sWi]);
state0.c_SO4   = 0.0000001*ones(G.cells.num,1);
state0.c_Ca    = 0.010000*ones(G.cells.num,1);
state0.c_Sr    = 0.002000*ones(G.cells.num,1);
state0.c_Ba    = 0.000100*ones(G.cells.num,1);
state0.c_Na    = 0.061200*ones(G.cells.num,1);
state0.c_Mg    = 0.001046*ones(G.cells.num,1);
state0.c_CO3   = 0.000427*ones(G.cells.num,1);
state0.c_Cl    = 0.180000*ones(G.cells.num,1);

model = OilWaterScaleModel(G, rock, fluid);
model.initial_props.poro = initial_phi;
model.initial_props.sWi = sWi;
model.initial_props.T = 93; %temperature
model.initial_props.P = initial_Pressure;
model.initial_props.perm = initial_k;
model.pvCorrector = ones(NX*NY*NZ,1);


simTime = 10*365*day;
nstep   = 800;
startsteps_dt = 1/4*day;
startsteps = repmat(startsteps_dt, 20, 1);
reststeps_dt = 2*day;
reststeps = repmat(reststeps_dt, nstep-20, 1);
timesteps = [startsteps;reststeps];


injRate = 100*stb/day;


W = verticalWell([],G,rock,ceil(NX/2),ceil(NY/2),1,...
                'Type', 'bhp', 'Val', 1500*psia, ...
                'Name', 'P1','comp_i',[0 1],'sign',-1,'radius',7*inch);       


% % % Injectors
W = verticalWell(W,G,rock,1,1,1,...
                  'Type', 'bhp', 'Val', 3500*psia,...
                  'Name','I1','comp_i',[1 0], 'sign',1, 'radius',7*inch);
W = verticalWell(W,G,rock,NX,1,1,...
                  'Type', 'bhp', 'Val', 3500*psia,...
                  'Name','I2','comp_i',[1 0], 'sign',1, 'radius',7*inch);
W = verticalWell(W,G,rock,1,NY,1,...
                  'Type', 'bhp', 'Val', 3500*psia,...
                  'Name','I3','comp_i',[1 0], 'sign',1, 'radius',7*inch);
W = verticalWell(W,G,rock,NX,NY,1,...
                  'Type', 'bhp', 'Val', 3500*psia,...
                  'Name','I4','comp_i',[1 0], 'sign', 1, 'radius',7*inch);
             
W(1).c_SO4=0;
W(2).c_SO4=0.01;
W(3).c_SO4=0.01;
W(4).c_SO4=0.01;
W(5).c_SO4=0.01;
W(1).c_Ca=0;
W(2).c_Ca=0.006;
W(3).c_Ca=0.006;
W(4).c_Ca=0.006;
W(5).c_Ca=0.006;
W(1).c_Sr=0;
W(2).c_Sr=0.000090;
W(3).c_Sr=0.000090;
W(4).c_Sr=0.000090;
W(5).c_Sr=0.000090;
W(1).c_Ba=0;
W(2).c_Ba=0;
W(3).c_Ba=0;
W(4).c_Ba=0;
W(5).c_Ba=0;
W(1).c_Na=0;
W(2).c_Na=0.013200;
W(3).c_Na=0.013200;
W(4).c_Na=0.013200;
W(5).c_Na=0.013200;
W(1).c_Mg=0;
W(2).c_Mg=0.001311;
W(3).c_Mg=0.001311;
W(4).c_Mg=0.001311;
W(5).c_Mg=0.001311;
W(1).c_CO3=0;
W(2).c_CO3=0.000156;
W(3).c_CO3=0.000156;
W(4).c_CO3=0.000156;
W(5).c_CO3=0.000156;
W(1).c_Cl=0;
W(2).c_Cl=0.021100;
W(3).c_Cl=0.021100;
W(4).c_Cl=0.021100;
W(5).c_Cl=0.021100;


schedule = simpleSchedule(timesteps,'W',W);

% [wellSols, states,rockProps, report] = simulateScheduleAD(state0, model, schedule);
fn = getPlotAfterStep(state0, model, schedule, ...
                      'plotReservoir', true, ...
                      'field', 'c_SO4',...
                      'wells',W,'view',[20 50]);
[wellSols, states,rockProps, report] = simulateScheduleAD(state0, model, schedule,...
                                    'afterStepFn', fn);

qOr   = deal(zeros(nstep,1));
qWr   = deal(zeros(nstep,1));
qw_SO4= deal(zeros(nstep,1));
qw_Ca = deal(zeros(nstep,1));
qw_Sr = deal(zeros(nstep,1));
qw_Ba = deal(zeros(nstep,1));
qw_Na = deal(zeros(nstep,1));
qw_Mg = deal(zeros(nstep,1));
qw_CO3= deal(zeros(nstep,1));
qw_Cl = deal(zeros(nstep,1));
qOt   = deal(zeros(nstep,1));
p_ave = deal(zeros(nstep,1));
phi   = deal(zeros(nstep,1));
perm_I1  = deal(cell(nstep,1));
PV = sum(G.cells.volumes)*initial_phi;
PV = convertTo(PV,stb);                                    
 
                                    
figure
grid
hold
for i=1:nstep
    qOr(i)        = wellSols{i,1}.qOs;
    qOr(i)        = abs(qOr(i));
    qWr(i)        = wellSols{i,1}.qWs;
    qWr(i)        = abs(qWr(i));
    [qw_SO4(i),~] = wellSols{i,1}.qW_SO4;
    qw_SO4(i)     = qw_SO4(i)*-1/qWr(i);
    [qw_Ca(i),~]  = wellSols{i,1}.qW_Ca;
    qw_Ca(i)      = qw_Ca(i) *-1/qWr(i);
    [qw_Sr(i),~]  = wellSols{i,1}.qW_Sr;
    qw_Sr(i)      = qw_Sr(i) *-1/qWr(i);
    [qw_Ba(i),~]  = wellSols{i,1}.qW_Ba;
    qw_Ba(i)      = qw_Ba(i) *-1/qWr(i);
    [qw_Na(i),~]  = wellSols{i,1}.qW_Na;
    qw_Na(i)      = qw_Na(i) *-1/qWr(i);
    [qw_Mg(i),~]  = wellSols{i,1}.qW_Mg;
    qw_Mg(i)      = qw_Mg(i) *-1/qWr(i);
    [qw_CO3(i),~]  = wellSols{i,1}.qW_CO3;
    qw_CO3(i)      = qw_CO3(i) *-1/qWr(i);
    [qw_Cl(i),~]  = wellSols{i,1}.qW_Cl;
    qw_Cl(i)      = qw_Cl(i) *-1/qWr(i);
end
t = 1:nstep;
plot(t,convertTo(qOr,(stb/day)),'*')
plot(t,convertTo(qWr,(stb/day)),'-')
xlabel('timestep')
ylabel('Rate (stb/day)')
legend('Oil Production Rate','Water Production Rate')  

figure
grid
hold
plot(t,qw_Ca*1e6,'b*')
plot(t,qw_Sr*1e6,'r-')
plot(t,qw_Ba*1e6,'y+')
plot(t,qw_SO4*1e6,'k.-')
plot(t,qw_Na*1e6,'gx')
plot(t,qw_Cl*1e6,'ms')
plot(t,qw_CO3*1e6,'b--')
plot(t,qw_Mg*1e6,'cd')
xlabel('time step')
ylabel('ppm')
legend('Ca','Sr','Ba','SO4','Na','Cl','CO3','Mg')

figure
for i=1:nstep
    p_ave(i) = mean(states{i,1}.pressure);
end
plot(t,convertTo(p_ave,psia))
xlabel('timestep')
ylabel('pressure (psia)')
grid

figure
for i=1:nstep
    phi(i)=mean(rockProps.poro{i,1});
end
plot(t,phi)
xlabel('timestep')
ylabel('phi')
grid

figure
for i=1:nstep
    a=rockProps.perm{i,1};
    perm_I1(i) = convertTo(a(1),milli*darcy);
end
plot(t,perm_I1)
xlabel('timestep')
ylabel('Injecting Grid Perm (mD)')
grid
