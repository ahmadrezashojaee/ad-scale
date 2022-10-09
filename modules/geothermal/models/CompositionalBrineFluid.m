classdef CompositionalBrineFluid < CompositionalMixture
    
    properties
        molecularDiffusivity;
    end
    
    methods
        %-----------------------------------------------------------------%
        function fluid = CompositionalBrineFluid(names, molarMass, molecularDiffusivity, varargin)
            fluid = fluid@CompositionalMixture(names, nan, nan, nan, nan, molarMass, varargin{:});
            fluid.molecularDiffusivity = molecularDiffusivity;
        end
    end
    
end