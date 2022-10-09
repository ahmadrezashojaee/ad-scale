classdef WaterModelDP < DualPorosityReservoirModel
    % Single phase water model. 
    properties

    end
    
    methods
        function model = WaterModelDP(G, rock, fluid, varargin)
            model = model@DualPorosityReservoirModel(G, rock, fluid);
            % Only enable water
            model.oil = false;
            model.gas = false;
            model.water = true;
            model.useCNVConvergence = false;
            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsWaterDP(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
    end
end

    %{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
