classdef DualContMechMechanicModel < DualContinuumMechanicModel
%
% SYNOPSIS:
%   model = DualContMechMechanicModel(G, rock, rock_matrix, mech_problem, varargin)
%
% DESCRIPTION: This model is for the mechanical part of a fully coupled
% model. It adds some few functionalities that are needed to couple the solver
% with a fluid model.
%
% PARAMETERS:
%   G            - Grid structure
%   rock         - Rock structure for fracture phase
%   rock_matrix  - Rock structure for matrix
%   mech_problem - Structure that contains the mechanical parameters of the system
%
% RETURNS:
%   class instance
%
% EXAMPLE: 
%
% SEE ALSO: DualContMechWaterModel  
%
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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

    properties
        % Primary variables that are used in the mechanic system
        primaryVarNames;
    end

    methods
        function model = DualContMechMechanicModel(G, rock, rock_matrix, mech_problem, varargin)
        % constructor
            model = model@DualContinuumMechanicModel(G, rock, rock_matrix, mech_problem, varargin{:});
            model.primaryVarNames = {'xd'};
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            error(['The DualContMechMechanicalModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from DualContMechFluidModel.'])
        end

        function  fds = getAllVarsNames(model) % getter function for the varsnames
        % list all the variables that are recognized and can be handled by the model
            fds = {'xd', 'uu', 'u', 'stress', 'strain', 'strain_frac', 'strain_mat', 'vdiv'};
        end

        function varnames = getAllPrimaryVariables(model) % getter function
        % list all the primary variables that are recognized and can be
        % handled by the model, used by the updateState member function.
            varnames = model.primaryVarNames;
        end

        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
            % updates the state variable that belongs to the model, that is
            % the mechanical variables. The fluid variables in the state will
            % be updated by the fluid model. 
            vars = problem.primaryVariables;
            ind = false(size(vars));
            mechvars = model.getAllPrimaryVariables();
            [lia, loc] = ismember(mechvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');

            ind(loc) = true;
            problem.primaryVariables = vars(ind);
            dx   = dx(ind);

            [state, report] = updateState@DualContinuumMechanicModel(model, state, problem, ...
                                                        dx, drivingForces);                        
            state = addDerivedQuantities(model, state);
        end

    end
end
