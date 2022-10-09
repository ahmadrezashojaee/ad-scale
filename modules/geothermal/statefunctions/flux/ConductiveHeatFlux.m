classdef ConductiveHeatFlux < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function hf = ConductiveHeatFlux(model, varargin)
            hf@StateFunction(model, varargin{:});
            hf = hf.dependsOn({'RockHeatTransmissibility', 'FluidHeatTransmissibility'});
            hf = hf.dependsOn({'Temperature'}, 'PVTPropertyFunctions');
            hf.label = 'H_c';
        end
        
        %-----------------------------------------------------------------%
        function H = evaluateOnDomain(prop,model, state)
            op = model.operators;
            [Tr, Tf] = prop.getEvaluatedDependencies(state, 'RockHeatTransmissibility' , ...
                                                            'FluidHeatTransmissibility');
            Tr = Tr(model.operators.internalConn);
            Tf = Tf(model.operators.internalConn);
            T = model.getProps(state, 'Temperature');
            H = -(Tr + Tf).*op.Grad(T);
        end 
    end
    
end

