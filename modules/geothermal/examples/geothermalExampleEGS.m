%% Enhanced geothermal system (EGS)
% This example shows how to simulate a small geothermal system in an
% artificial fracture newtork. The grid is constructed by extruding a 2D
% PEBI grid with refinement around volumetric fractures vertically.

%% Add necessary MRST modules
mrstModule add ad-core ad-props ad-blackoil geothermal compositional upr ...
    mrst-gui example-suite
mrstVerbose on

%% Set up example
example = MRSTExample('small_egs_geothermal');

%% Plot setup
example.figure();
plotGrid(example.model.G, example.model.G.cells.tag, 'faceColor', [1,1,1]*0.8, 'edgeColor', 'none');
plotGrid(example.model.G, 'faceColor', 'none', 'edgeAlpha', 0.1);
example.setAxisProperties(gca);
camlight()
plotWell(example.model.G, example.schedule.control(1).W, 'color', 'k', 'fontSize', 30);
axis off

%% Simulate
problem = example.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem, 'prompt', true);
simulatePackedProblem(problem);

%% Interactive plot of results
close all
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
example.plot(states);
colormap(hot);
plotWellSols(wellSols, example.schedule.step.val)

%% Temperature at selected timesteps
T = getWellOutput(wellSols, 'T');
Tmin = min(min(T));
Tmax = max(max(T));
steps = [8, 15, 30];
a = 0.6; cmap = hot.*a + (1-a);
for i = 1:numel(steps)
   example.figure();
   plotCellData(example.model.G, states{steps(i)}.T, 'edgeAlpha', 0.2);
   example.setAxisProperties(gca);
   plotWell(example.model.G, example.schedule.control(1).W, 'color', 'k', 'fontSize', 0.01);
   colormap(cmap), camlight();
   axis off;
   caxis([Tmin, Tmax]);
end

%% Plot EGS efficiency
close all
p   = getWellOutput(wellSols, 'bhp');
T   = getWellOutput(wellSols, 'T');
q   = getWellOutput(wellSols, 'qWs');

[h, rho] = deal(zeros(size(p)));
for i = 1:2
    h(:, i)   = example.model.fluid.hW(p(:,i), T(:,i));
    rho(:, i) = example.model.fluid.rhoW(p(:,i), T(:,i));
end
qH  = abs(q.*rho.*h);

eff = (qH(:,2))./qH(:,1);
time = cumsum(example.schedule.step.val);

figure('Position', [0, 0, 800, 200])

hold on
kWatt = kilo*joule/second;
plot(time/year, eff, 'color', 'k', 'linew', 2);
% Indicate timesteps plotted above by circles
plot(time(steps)/year, eff(steps), 'ok', 'linew', 2);
axis([[time(5), time(end)]/year, min(eff(5:end))*0.95, max(eff(5:end))*1.05]);
set(gca, 'Box', true, 'FontSize', 13);
xlabel('Time (years)')