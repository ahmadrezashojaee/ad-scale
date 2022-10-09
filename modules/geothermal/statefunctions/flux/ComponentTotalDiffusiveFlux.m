classdef ComponentTotalDiffusiveFlux < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function df = ComponentTotalDiffusiveFlux(model, varargin)
            df@StateFunction(model, varargin{:});
            df = df.dependsOn({'MolecularTransmissibility'});
            df = df.dependsOn({'ComponentPhaseMassFractions', 'Density'}, 'PVTPropertyFunctions');
            df.label = 'D_i';
        end
        
        %-----------------------------------------------------------------%
        function D = evaluateOnDomain(prop, model, state)
            op = model.operators;
            T = prop.getEvaluatedDependencies(state, 'MolecularTransmissibility');
            ix       = ~cellfun(@isempty, T);
            T(ix)    = cellfun(@(T) T(model.operators.internalConn), T(ix), 'UniformOutput', false);
            [X, rho] = model.getProps(state, 'ComponentPhaseMassFractions', ...
                                             'Density'                    );
            ncomp = model.getNumberOfComponents();
            nph   = model.getNumberOfPhases();
            D     = cell(ncomp,1);
            [D{:}] = deal(0);
            for ph = 1:nph
                rhof = model.operators.faceAvg(rho{ph});
                for c = 1:ncomp
                    if isempty(T{c})
                        continue
                    end
                    D{c} = D{c} - rhof.*T{c}.*op.Grad(X{c, ph});
                end
            end
        end 
    end
    
end