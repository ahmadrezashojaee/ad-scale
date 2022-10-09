%% Subdomain example
% In this example, we show how to construct a submodel of a full model and
% solve the corresponding subproblem
mrstModule add ad-core ad-props ad-blackoil sequential spe10 ...
    mrst-gui example-suite

%% Get example
% Get SPE10 subset
example = MRSTExample('spe10_wo', 'layers', 38, 'J', 1:110);
% Replace wells
W0 = example.schedule.control(1).W;
well = @(W, i, j, varargin) verticalWell(W, example.model.G, example.model.rock, i, j, [], ...
                      'Radius', W0(end).r, 'WI', W0(end).WI, 'comp_i', [1,0], varargin{:});
W = well([], 51,   1, 'type', 'rate', 'val', W0(end).val);
W = well(W , 32, 110, 'type', 'bhp' , 'val', W0(1).val  );
% Update example
example.schedule.control(1).W = W;
example.schedule.step.val = example.schedule.step.val(1:40);
example.schedule.step.control = example.schedule.step.control(1:40);
example.name = 'spe10-layer-38';
example.axisProperties.PlotBoxAspectRatio = [1,1,1];
example.figureProperties.Size = [600,600];
% Plot setup
example.plot(example.model.rock, 'log10', true);

%% Set up transport problem
% We compute a representative flow field, and solve the transport
% subproblem only
modelSeq = getSequentialModelFromFI(example.model);
% Upstate initial state with flow properties
state0 = standaloneSolveAD(example.state0, modelSeq.pressureModel, 100*day, 'W', W);
state0 = rmfield(state0, 'statePressure');
state0 = rmfield(state0, 'FlowProps');
% Update example
example.model = modelSeq.transportModel;
example.state0 = state0;

%% Simulate transport
problem = example.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem, 'prompt', true);
simulatePackedProblem(problem);

%% Inspect result
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
example.plot(states);

%% Identify high-flow regions
flux = states{end}.flux;
v    = sqrt(sum(faceFlux2cellVelocity(example.model.G, sum(flux,2)).^2,2));
example.plot(v, 'log10', true)
cells = log10(v) > -7.5;
% Add padding of three cells
C = getConnectivityMatrix(example.model.parentModel.operators.N);
for i = 1:3
    cells = cells | C*cells;
end
plotGrid(example.model.G, cells, 'faceColor', 'none')

%% Construct a submodel of the high-flow region only
mrstModule add domain-decomposition
submodel = SubdomainModel(example.model, cells);
% Make subexample
subexample          = example;
subexample.model    = submodel;
subexample.state0   = getSubState(subexample.state0, submodel.mappings);
subexample.schedule = getSubSchedule(subexample.schedule, submodel.mappings);
subexample.name     = [subexample.name, '-sub'];
% The submodel grid structure is even plottable by itself
subexample.plot(subexample.model.parentModel.parentModel.rock, 'log10', true);

%% Simulate subproblem
subproblem = subexample.getPackedSimulationProblem();
clearPackedSimulatorOutput(subproblem, 'prompt', true);
simulatePackedProblem(subproblem);

%% Inspect results
[wellSolsSub, statesSub, reportsSub] = getPackedSimulatorOutput(subproblem);
subexample.plot(statesSub);

%% Plot difference
% Map states to full state
map = submodel.mappings; map.cells.fields = {'s'};
statesSub2Full = cellfun(@(s) mapState(example.state0, s, map, 'mapWellSol', false), ...
                                       statesSub, 'uniformOutput', false);
getsw = @(state, c) state.s(c,1);
figure('Position', [0, 0, 1000, 500])
cmap = flipud(winter); cmap(1,:) = [1,1,1];
[bf, c] = boundaryFaces(example.model.G);
vf = example.model.G.faces.normals(bf,3)./example.model.G.faces.areas(bf) < 1e-10;
% Plot saturation
subplot(1,2,1), title('Reference water saturation');
h = plotCellData(example.model.G, getsw(states{1}, ':'), 'edgeColor', 'none');
example.setAxisProperties(gca); caxis([0.208, 0.8]); colorbar('Location', 'southoutside');
plotGrid(subexample.model.G, 'faceColor', 'none', 'edgeAlpha', 0.2); plotFaces(example.model.G, bf(vf));
% Plot difference
subplot(1,2,2), title('Difference');
hd = plotCellData(example.model.G, abs(getsw(states{1}, ':') - getsw(statesSub2Full{1}, ':')), 'edgeColor', 'none');
example.setAxisProperties(gca); caxis([0, 0.201]); colorbar('Location', 'southoutside');
plotGrid(subexample.model.G, 'faceColor', 'none', 'edgeAlpha', 0.2); plotFaces(example.model.G, bf(vf));
colormap(cmap);
% Loop through timesteps
for i = 1:numel(statesSub2Full)
    h.FaceVertexCData = getsw(states{i}, c);
    hd.FaceVertexCData = abs(getsw(states{i}, c) - getsw(statesSub2Full{i}, c));
    pause(0.2);
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
