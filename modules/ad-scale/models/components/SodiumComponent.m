classdef SodiumComponent < ConcentrationComponent
    properties

    end

    methods
        function c = SodiumComponent()
            % sulfate transport in water phase
            wIx = 1;
            c@ConcentrationComponent('sodium', wIx);
            
            c = c.functionDependsOn('getComponentDensity', 'sodium',        'state');
            c = c.functionDependsOn('getComponentDensity', 'ShrinkageFactors',  'PVTPropertyFunctions');

            c = c.functionDependsOn('getComponentMass', {'s', 'sodium'},    'state');
            c = c.functionDependsOn('getComponentMass', 'ShrinkageFactors',     'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMass', 'PoreVolume',           'PVTPropertyFunctions');
            
            c = c.functionDependsOn('getComponentMobility', 'sodium',       'state');
            c = c.functionDependsOn('getComponentMobility', 'ShrinkageFactors', 'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMobility', 'Mobility',         'FlowPropertyFunctions');
        end

        function c = getComponentDensity(component, model, state, varargin)
            c_Na = model.getProp(state, 'sodium');
            b = model.getProps(state, 'ShrinkageFactors');
            nph = numel(b);
            c = cell(1, nph);
            % water phase index
            wIx = 1;
            c{wIx} = c_Na .* b{wIx};
        end


        function c = getComponentMass(component, model, state, varargin)
             f = model.fluid;
             c_Na = model.getProp(state, 'sodium');
             pv = model.getProp(state, 'PoreVolume');
             b = model.getProps(state, 'ShrinkageFactors');

             nph = model.getNumberOfPhases;
             c = cell(1, nph);

             wIx = 1;
             bW = b{wIx};
             sw = model.getProp(state, 'sW');
             % In mobile water
             acc = sw.*c_Na.*bW;             

             c{wIx} = pv.*acc;
        end

   
        function cmob = getComponentMobility(component, model, state, varargin)

             mob = model.getProps(state, 'Mobility');
             wIx = 1;
             mobW = mob{wIx};

             nphase = model.getNumberOfPhases;
             cmob = cell(1, nphase);
             cmob{wIx} = mobW;
        end
        
        
        function c = getInjectionMassFraction(component, model, force)
            c = vertcat(force.c_Na)./model.fluid.rhoWS;
        end

    end
end

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
