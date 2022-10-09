classdef FluidThermalConductivity < StateFunction
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = FluidThermalConductivity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'Temperature', 'PoreVolume'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn('pressure', 'state');
            gp.label = '\Lambda_F';
        end
        
        %-----------------------------------------------------------------%
        function lambdaF = evaluateOnDomain(prop, model, state)
            [p, T, pv] = model.getProps(state, 'pressure', 'Temperature', 'PoreVolume');
            lambdaF = prop.evaluateFluid(model, 'lambdaF', p, T);
            vol = model.G.cells.volumes;
            lambdaF = lambdaF.*pv./vol;
        end
    end
    
end
