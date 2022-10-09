classdef ComponentPhaseMoleFractionsBrine < StateFunction

        properties
        end
        
        methods
            %-------------------------------------------------------------%
            function gp = ComponentPhaseMoleFractionsBrine(model, varargin)
                gp@StateFunction(model, varargin{:});
                gp = gp.dependsOn({'x'}, 'state');
                gp.label = 'x_{i,\alpha}';
            end
            
            %-------------------------------------------------------------%
            function x = evaluateOnDomain(prop, model, state) %#ok
                x = model.getProps(state, 'x');
                if ~iscell(x)
                    x = expandMatrixToCell(x);
                end
                x = reshape(x, [], 1);
            end
        end
        
end