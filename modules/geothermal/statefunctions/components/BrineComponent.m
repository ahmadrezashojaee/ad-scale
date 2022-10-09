classdef BrineComponent < GenericComponent
    
    properties
        componentIndex
        molecularDiffusivity
        surfacePhaseMassFractions = 1;
    end
    
    methods
        %-----------------------------------------------------------------%
        function component = BrineComponent(name, molarMass, diffusivity, index)
            component                      = component@GenericComponent(name);
            component.molarMass            = molarMass;
            component.molecularDiffusivity = diffusivity;
            component.componentIndex       = index;
        end
        
        %-----------------------------------------------------------------%
        function componetDensity = getComponentDensity(component, model, state)
            % Density of component in each phase (mass per unit of volume)
            [rho, X] = model.getProps(state, 'Density', 'ComponentPhaseMassFractions');
            ix = strcmp(model.getComponentNames, component.name);
            componetDensity = {rho{1}.*X{ix}};
        end
        
        %-----------------------------------------------------------------%
        function c = getPhaseComponentFractionInjection(component, model, state, force)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            comp_i = vertcat(force.components);
            comp_i = model.getMassFraction(comp_i);
            index = component.componentIndex;
            ci = comp_i(:, index);
            if any(ci ~= 0)
                for i = 1:nph
                    c{i} = ci;
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            for i = 1:nph
                c{i} = component.surfacePhaseMassFractions(i);
            end
        end
    end
    
end
