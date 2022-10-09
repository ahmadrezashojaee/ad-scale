function problem = getSubSchedule(problem, steps)
% Example function for optional QoI class property `processProblemFn`. This
% function extracts a subset of the schedule after setting the sample.
% Useful for e.g. when the schedule read from a deck.

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    schedule = problem.SimulatorSetup.schedule;
    schedule.step.val     = schedule.step.val(steps);
    schedule.step.control = schedule.step.control(steps);
    problem.SimulatorSetup.schedule = schedule;
end
