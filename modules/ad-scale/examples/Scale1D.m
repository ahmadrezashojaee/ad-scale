mrstModule add ad-core ad-blackoil ad-scale ad-props

% Define grid
[NX,NY,NZ]=deal(20,1,1);
[dx,dy,dz]=deal(1000,10,5);
G=computeGeometry(cartGrid([NX,NY,NZ],[dx,dy,dz]));
initial_k = 100*milli*darcy;
initial_phi = 0.2;

rock = makeRock(G, initial_k, initial_phi);

ref_depth=1000;
fluid = initSimpleADIFluid('mu',    [1.2, 5]*centi*poise, ...
                           'rho',   [1019, 750]*kilogram/meter^3, ...
                           'n',     [2, 2], ...
                           'cR',    10e-6/psia, ...
                           'c',     [1e-6*psia^-1 9e-6*psia^-1],...
                           'smin',  [0.2 0.2],...
                           'phases', 'wo');
                       
% Construct reservoir model
gravity reset on
sWi = 0.2;
initial_Pressure = 3000*psia;
state0      = initResSol(G, initial_Pressure, [sWi, 1-sWi]);
state0.c_SO4    = 0.000001*ones(G.cells.num,1);
state0.c_Ca    = 0.006*ones(G.cells.num,1);
state0.c_Sr    = 0.002000*ones(G.cells.num,1);
state0.c_Ba    = 0.000100*ones(G.cells.num,1);
state0.c_Na    = 0.061200*ones(G.cells.num,1);
state0.c_Mg    = 0.001046*ones(G.cells.num,1);
state0.c_CO3    = 0.000427*ones(G.cells.num,1);
state0.c_Cl    = 0.18*ones(G.cells.num,1);

model = OilWaterScaleModel(G, rock, fluid);
model.initial_props.poro = initial_phi;
model.initial_props.sWi = sWi;
model.initial_props.T = 90; %temperature
model.pvCorrector = ones(NX*NY*NZ,1);

%schedule
simTime = 5*365*day;
nstep   = 1000;
startsteps_dt = 3*day;
startsteps = repmat(startsteps_dt, 240, 1);
reststeps_dt = 5*day;
reststeps = repmat(reststeps_dt, nstep-240, 1);
timesteps = [startsteps;reststeps];



pv      = poreVolume(G, rock);
injRate = 100*stb/day;


W = verticalWell([],G,rock,NX,NY,1:NZ,...
                'Type', 'bhp', 'Val', 1500*psia, ...
                'Name', 'P1','comp_i',[0 1],'sign',-1, 'radius', 0.2);       


% % % Injectors
W = verticalWell(W,G,rock,1,1,1:NZ,...
                  'Type', 'bhp', 'val', 3500*psia,...
                  'Name','I1','comp_i',[1 0],'sign', 1, 'radius', 0.2);
               
W(1).c_SO4=0;
W(2).c_SO4=0.01;
W(1).c_Ca=0;
W(2).c_Ca=0;
W(1).c_Sr=0;
W(2).c_Sr=0;
W(1).c_Ba=0;
W(2).c_Ba=0;
W(1).c_Na=0;
W(2).c_Na=0.006000;
W(1).c_Mg=0;
W(2).c_Mg=0;
W(1).c_CO3=0;
W(2).c_CO3=0;
W(1).c_Cl=0;
W(2).c_Cl=0;

schedule = simpleSchedule(timesteps,'W',W);%, 'bc', bc);%,'W', W);
fn = getPlotAfterStep(state0, model, schedule, ...
                      'plotReservoir', true, 'view', [20, 8], ...
                      'field', 'c_Ca', 'Wells',W);
[wellsol, states, rockProps,report] = simulateScheduleAD(state0, model, schedule, ...
                                'afterStepFn', fn);
                                 
%% 
%plot oil rate
%plot oil rate
qOr   = deal(zeros(nstep,1));
qWr   = deal(zeros(nstep,1));
ppm_SO4 = deal(zeros(nstep,1));
ppm_Ca = deal(zeros(nstep,1));
ppm_Sr = deal(zeros(nstep,1));
ppm_Ba = deal(zeros(nstep,1));
ppm_Na = deal(zeros(nstep,1));
ppm_Cl = deal(zeros(nstep,1));
ppm_Mg = deal(zeros(nstep,1));
ppm_CO3 = deal(zeros(nstep,1));
p_ave = deal(zeros(nstep,1));
phi   = deal(zeros(nstep,1));
perm  = deal(cell(nstep,1));
Qo_T  = deal(zeros(nstep,1));
OIP   = convertTo(dx*dy*dz*initial_phi*(1-sWi),centi^3);
figure
for i=1:nstep
    qOr(i)   = states{i,1}.flux(NX+1,2);
    qOr(i)   = convertTo(qOr(i),centi^3/hour);
    qWr(i)   = states{i,1}.flux(NX+1,1);
    qWr(i)   = convertTo(qWr(i),centi^3/hour);
    ppm_SO4(i) = states{i,1}.c_SO4(NX)*1e6;
    ppm_Ca(i)  = states{i,1}.c_Ca(NX)*1e6;
    ppm_Sr(i)  = states{i,1}.c_Sr(NX)*1e6;
    ppm_Ba(i)  = states{i,1}.c_Ba(NX)*1e6;
    ppm_Na(i)  = states{i,1}.c_Na(NX)*1e6;
    ppm_Mg(i)  = states{i,1}.c_Mg(NX)*1e6;
    ppm_CO3(i) = states{i,1}.c_CO3(NX)*1e6;
    ppm_Cl(i)  = states{i,1}.c_Cl(NX)*1e6;
    

end
t = 1:nstep;
% t = t.*dt./3600;
plot(t,qOr)
xlabel('time step')
ylabel('Oil Rate (cc/hour)')
grid on
title('scale')

figure
hold on
grid on
plot(t,ppm_Ca,'*')
plot(t,ppm_SO4,'.-')
plot(t,ppm_Sr,'+')
plot(t,ppm_Ba,'-')
plot(t,ppm_Na,'o')
plot(t,ppm_Mg,'h')
plot(t,ppm_Cl,'d')
plot(t,ppm_CO3,'--')

ylabel('ppm')
xlabel('time step')
legend('Ca','SO4','Sr','Ba','Na','Mg','Cl','CO3')

figure
plot(t,100.*Qo_T./OIP,'-.b')
xlabel('time step')
ylabel('Rf (%)')
grid on
title('scale')


figure
for i=1:nstep
    p_ave(i) = mean(states{i,1}.pressure);
end
plot(t,convertTo(p_ave,psia))
xlabel('timestep')
ylabel('pressure (psia)')
grid
title('scale')


% figure
% plotToolbar(G,rockProps.poro)

figure
for i=1:nstep
    perm{i,1}=rockProps.perm{i,1};
    perm{i,1}=convertTo(perm{i,1},milli*darcy);
end
plotToolbar(G,perm)
xlabel('timestep')
ylabel('k')
grid on
title('scale')

% figure
% plotToolbar(G,rockProps.poro)
% title('scale')
%%