classdef AdvectiveHeatFlux < StateFunction & UpwindProperty
    
    properties (Access = protected)
        upwind_name; % Name of state function where upwind flag comes from
    end
    
    methods
        %-----------------------------------------------------------------%
        function hf = AdvectiveHeatFlux(model, varargin)
            if nargin < 2
                upstr = UpwindFunctionWrapperDiscretization(model);
            end
            if nargin < 3
                upwind_name = 'PhaseUpwindFlag';
            end
            hf@StateFunction(model);
            hf@UpwindProperty(upstr)
            hf.upwind_name = upwind_name;
            hf = hf.dependsOn({upwind_name, 'PhaseFlux'});
            hf = hf.dependsOn({'Enthalpy', 'Density'}, 'PVTPropertyFunctions');
            hf.label = 'H_{a,\alpha}';
        end
        
        %-----------------------------------------------------------------%
        function H = evaluateOnDomain(prop, model, state)
            [v, flag] = prop.getEvaluatedDependencies(state, 'PhaseFlux', prop.upwind_name);
            [h, rho]  = model.getProps(state, 'Enthalpy', 'Density');
            nph = numel(rho);
            H = cell(1, nph);
            for i = 1:nph
                H{i} = prop.faceUpstream(model, state, flag{i}, rho{i}.*h{i}).*v{i};
            end
        end 
    end
    
end

