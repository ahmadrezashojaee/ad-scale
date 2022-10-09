function src = getHeatFluxBoundary(model, state, src, drivingForces)
% Add heat flux from boundary conditions
    % Early return if no BCs are given
    if isempty(src.bc.sourceCells)
        src.bc.heatFlux = [];
        return
    end
    bc = drivingForces.bc;
    % Get reservoir properties
    [hr, pr, Tr, rhor, mobr, Tf, ThR, ThF] = getReservoirProperties(model, state, bc);
    % Get boundary face properties
    [Tbc, Xbc] = getBoundaryProperties(model, Tr, ThR, ThF, bc);
    qAdv  = computeAdvectiveHeatFlux(model, mobr, rhor, pr, hr, Tbc, Xbc, Tf, bc, src);
    qCond = computeConductiveHeatFlux(ThR, ThF, Tr, Tbc, bc);
    src.bc.heatFlux = qAdv + qCond;
end

%-------------------------------------------------------------------------%
function q = computeAdvectiveHeatFlux(model, mobr, rhor, pr, hr, Tbc, Xbc, Tf, bc, src)
    nph = model.getNumberOfPhases();
    q = 0;
    for i = 1:nph
        inflow = src.bc.phaseMass{i} > 0;
        pbc = getBoundaryPressure(model, mobr{i}, rhor{i}, Tf, pr, bc);
        hbc = model.fluid.hW(pbc, Tbc, Xbc);
        h = inflow.*hbc + ~inflow.*hr{i};
        q = q + src.bc.phaseMass{i}.*h;
    end
end

%-------------------------------------------------------------------------%
function q = computeConductiveHeatFlux(ThR, ThF, Tr, Tbc, bc)
    is_Hflux = ~isnan(bc.Hflux);
    q = -(ThR + ThF).*(Tr - Tbc);
    q(is_Hflux) = bc.Hflux(is_Hflux);
end

%-------------------------------------------------------------------------%
function [h, p, T, rho, mob, Tf, ThR, ThF]  = getReservoirProperties(model, state, bc)
    [h, p, T, rho, mob, ThR, ThF] = model.getProps(state, 'Enthalpy', ...
                                                          'pressure', ...
                                                          'Temperature', ...
                                                          'Density', ...
                                                          'Mobility', ...
                                                          'FluidHeatTransmissibility', ...
                                                          'RockHeatTransmissibility');
    faces = bc.face;
    cells = sum(model.G.faces.neighbors(faces,:),2);
    % Extract BC cell values
    getBCVal = @(v) cellfun(@(v) v(cells), v, 'UniformOutput', false);
    p   = p(cells);
    T   = T(cells);
    h   = getBCVal(h);
    rho = getBCVal(rho);
    mob = getBCVal(mob);
    % Extract BC face values
    Tf  = model.operators.T_all(faces);
    ThR = ThR(faces);
    ThF = ThF(faces);
end

%-------------------------------------------------------------------------%
function [Tbc, Xbc] = getBoundaryProperties(model, T, ThR, ThF, bc)
    Tbc = bc.T;
    cnames = model.getComponentNames();
    ix = strcmpi(cnames, 'NaCl');
    Xbc = model.getMassFraction(bc.components);
    if any(ix), Xbc = Xbc(:,ix); else, Xbc = []; end 
    is_Hflux = ~isnan(bc.Hflux);
    if any(is_Hflux)
        Tbc = model.AutoDiffBackend.convertToAD(Tbc, T);
        dT            = bc.Hflux(is_Hflux)./(ThR(is_Hflux) + ThF(is_Hflux));
        Tbc(is_Hflux) = T(is_Hflux) + dT;
    end
end

%-------------------------------------------------------------------------%
function p = getBoundaryPressure(model, mob, rho, T, p, bc)
    is_pressure = reshape(strcmpi(bc.type, 'pressure'), [], 1);
    is_flux     = reshape(strcmpi(bc.type, 'flux'), [], 1);
    % Get pressure at boundaries open to flow
    if any(is_pressure)
        p(is_pressure) = bc.value(is_pressure);
    end
    cells = sum(model.G.faces.neighbors(bc.face,:), 2);
    if any(is_flux)
        % Flux BCs requires reconstruction. Assuming simple
        G = model.G;
        if any(strcmpi(G.type,'topSurfaceGrid'))
            dzbc = model.gravity(3)*(G.cells.z(cells) - G.faces.z(bc.face));
        else
            g    = model.getGravityVector();
            dz   = G.cells.centroids(cells,:) - G.faces.centroids(bc.face,:);
            dzbc = dz*g';
        end
        dp         = bc.value(is_flux)./(mob(is_flux).*T(is_flux));
        p(is_flux) = p(is_flux) + dp - rho(is_flux).*dzbc(is_flux);
    end
end