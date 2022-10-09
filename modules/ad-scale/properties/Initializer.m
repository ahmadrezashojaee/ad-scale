function state0 = Initializer(model,state0)
dt = 10000*day;
nstep = 1;
timesteps = repmat(dt, nstep, 1);
schedule = simpleSchedule(timesteps);
[~,states,~,~] = simulateScheduleAD(state0,model,schedule);
clc
state0=states{nstep,1};
end
