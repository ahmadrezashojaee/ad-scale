classdef DynamicTransmissibility < StateFunction
    % Transmissibility for internal faces. May include an optional
    % pressure-dependent multiplier from a field in the fluid model.
    properties
        harmonicAvgOperator
        twoPointOperator
        conductivity_name
    end
    
    methods
        %-----------------------------------------------------------------%
        function dt = DynamicTransmissibility(model, conductivity_name)
            if nargin < 2
                conductivity_name = 'Permeability';
            end
            dt@StateFunction(model);
            dt = dt.dependsOn(conductivity_name);
            dt.conductivity_name   = conductivity_name;
            dt.twoPointOperator    = getTwoPointOperator(model.G);
            dt.harmonicAvgOperator = getHarmonicAvgOpeartor(model.G);
            dt.label = 'T';
        end
        
        %-----------------------------------------------------------------%
        function T = evaluateOnDomain(prop, model, state)
            lambda = prop.getEvaluatedDependencies(state, prop.conductivity_name);
            if iscell(lambda)
                T = cell(numel(lambda),1);
                for i = 1:numel(lambda)
                    if any(lambda{i})
                        T{i} = prop.getTransmissibility(lambda{i});
                    end
                end
            else
                T = prop.getTransmissibility(lambda);
            end
        end
        
        %-----------------------------------------------------------------%
        function T = getTransmissibility(prop, lambda)
            T = prop.twoPointOperator(lambda);
            T = prop.harmonicAvgOperator(T);
        end
    end
end

%-------------------------------------------------------------------------%
function tp = getTwoPointOperator(G)
    % Mappings from cells to its faces
    cells = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    faces = G.cells.faces(:,1);
    % Vector from cell to face centroid
    C = G.faces.centroids(faces,:) - G.cells.centroids(cells,:);
    % Oriented normals
    sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1;
    N   = bsxfun(@times, sgn, G.faces.normals(faces, :));
    % Make function
    cn  = sum(C.*N,2)./sum(C.*C,2);
    tp = @(lambda) cn.*lambda(cells);
end

%-------------------------------------------------------------------------%
function ha = getHarmonicAvgOpeartor(G)
    % Harmonig averaging operator
    faces = G.cells.faces(:,1);
    M = sparse(faces, 1:numel(faces), 1, G.faces.num, numel(faces));
    ha = @(T) 1./(M*(1./T));
end