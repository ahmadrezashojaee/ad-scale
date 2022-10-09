%% Introduction to nonlinear domain decomposition in MRST
% This example gives an introduction to the domain decomposition module,
% and shows how it can be used to solve a problem using additive nonlinear
% domain decomposition
mrstModule add ad-core ad-props ad-blackoil coarsegrid example-suite
mrstVerbose on

%% Get example
% We consider the ''qfw_peaks_wo'' example from the example-suite module.
% This is an inverted five-spot pattern on a permeability field made using
% tiles of the built-in MATLAB function peaks, creating a beautiful
% displacement front pattern
mrstModule add example-suite mrst-gui
n = 47;
example = MRSTExample('ifs_peaks_wo', 'ncells', n, 'tiles', 2);
example.plot(example.model.rock, 'log10', true);

%% Submodels
% Solving one timestep with additive nonlinear domain decomposition amounts
% to decompose the full domain into a number of subdomains, and solve the
% full timestep in each subdomain with unknowns in all other subdomains
% fixed. This is repeated until all subdomains are converged
% simultaneously. We construct a subdomain problem by means of a
% SubdomainModel that wrapps around the model. We illustrate this by making
% a subdomain model consisting of the third quadrant of the full model.
mrstModule add domain-decomposition
[ii, jj] = gridLogicalIndices(example.model.G);
cells = ii > n & jj > n;
submodel = SubdomainModel(example.model, cells);
% This model has a field 'mappings' that contains all necessary mappings
% for cells and faces.
colors = lines(2);
map    = submodel.mappings;
pargs  = {'edgeAlpha', 0.1};
example.figure();
plotGrid(example.model.G, 'faceColor', 'none', pargs{:});
h= [];
h(1) = plotGrid(example.model.G, map.cells.internal, 'faceColor', colors(1,:), pargs{:});
h(2) = plotGrid(example.model.G, map.cells.external, 'faceColor', colors(2,:), pargs{:});
example.setAxisProperties(gca);
legend(h, {'Internal', 'external'}, 'Location', 'northwest');

%% Get subdomain problem
% For convenience, we make a subexample as well. We use getSubState and
% getSubSchedule to get the corresponding submodel initial state and
% schedule.
subexample          = example;
subexample.name     = [subexample.name, '-submodel'];
subexample.model    = submodel;
subexample.state0   = getSubState(subexample.state0, submodel.mappings);
subexample.schedule = getSubSchedule(subexample.schedule, submodel.mappings);
% The submodel grid structure is even plottable by itself
subexample.plot(subexample.model.parentModel.rock, 'log10', true);

%% Simulate subdomain
% For illustrational purposes, we impose no-flow boundary conditions on the
% subdomain interfaces that are internal
subexample.model.noflowBC = true;
subproblem = subexample.getPackedSimulationProblem();
clearPackedSimulatorOutput(subproblem, 'prompt', true);
simulatePackedProblem(subproblem);

%% Inspect resutls
[~, substates] = getPackedSimulatorOutput(subproblem);
subexample.plot(substates);

%% Partition domain
% In order to do domain decomposition, we partition the domain using
% partitionCartGrid from the coarsegrid module to construct a rectilinear
% 5x5 subdomain partition
mrstModule add coarsegrid
partition = partitionCartGrid(example.model.G.cartDims, [5,5]);
example.plot(partition); % Plot subdomain partition

%% Set up nonlinear domain decomposition model
% Additive NLDD tends to converge poorly for elliptic problems. We therfore
% use sequential splitting, and use additive NLDD for the transport
% subproblem.
mrstModule add sequential
modelSeq = getSequentialModelFromFI(example.model);
% Make sequential example
exampleSeq       = example;
exampleSeq.model = modelSeq;
exampleSeq.name  = [exampleSeq.name, '-seq'];
% Nonlinear domain decompositioning is implemented in the
% DomainDecompositionModel. This takes the parent model and partition as
% input arguments. Each subdomain can be solved concurrently, so we do this
% if we have the parallel computing toolbox at hand.
parallel = ~isempty(ver('parallel'));
modelSeq.transportModel = DomainDecompositionModel(modelSeq.transportModel, partition, 'parallel', parallel);
% Make example
exampleSeqDD       = exampleSeq;
exampleSeqDD.model = modelSeq;
exampleSeqDD.name  = [exampleSeqDD.name, '-nldd'];

%% Simulate with NLDD in transport
problemSeqDD = exampleSeqDD.getPackedSimulationProblem();
clearPackedSimulatorOutput(problemSeqDD, 'prompt', true);
simulatePackedProblem(problemSeqDD);

%% Simulate withouth NLDD for reference
problemSeq = exampleSeq.getPackedSimulationProblem();
clearPackedSimulatorOutput(problemSeq, 'prompt', true);
simulatePackedProblem(problemSeq);

%% Inspect results
% Hint: The number of iterations used in each subdomain are stored in
% state.iterations
[wellSolsSeqDD, statesSeqDD, reportsSeqDD] = getPackedSimulatorOutput(problemSeqDD);
[wellSolsSeq  , statesSeq  , reportsSeq  ] = getPackedSimulatorOutput(problemSeq);
exampleSeqDD.plot(statesSeqDD);
plotWellSols({wellSolsSeq, wellSolsSeqDD});

%% Compare iterations
% Resolving nonlinearities locally can be very beneficial. We compare the
% number of nonlinear iterations used to solve the transport subproblem for
% each solver
its   = getReportOutput(reportsSeq  , 'solver', 'TransportSolver');
itsDD = getReportOutput(reportsSeqDD, 'solver', 'TransportSolver');
figure('Position', [0, 0, 800, 400]); hold on
plot(cumsum(its.total)  , 'LineWidth', 2)
plot(cumsum(itsDD.total), 'LineWidth', 2)
hold off, box on;
xlim([1, numel(its.total)]);
xlabel('Timestep')
title('Cumulative transport iterations')
legend({'Global', 'Domain decomposition'}, 'location', 'north west')

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
