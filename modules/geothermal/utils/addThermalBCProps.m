function bc = addThermalBCProps(bc, varargin)
% Add thermal boundary conditions to existing bc struct. Can be either a
% fixed temperature, or a heat flux.
    opt = struct('T'    , nan, ...
                 'Hflux', nan);
    opt = merge_options(opt, varargin{:});
    opt = checkInput(bc, opt, 'T'    );
    opt = checkInput(bc, opt, 'Hflux');
    assert(~any(opt.T > 0 & ~isnan(opt.Hflux)), ...
                           'Multiple thermal BCs given for the same face');
    assert(~any(isnan(opt.T) & isnan(opt.Hflux)), ...
           ['Please provide either a temperature or heat flux for ', ...
            'each boundary face']                                  );
    bc.T     = opt.T;
    bc.Hflux = opt.Hflux;
end

%-------------------------------------------------------------------------%
function opt = checkInput(bc, opt, name)
    nf = numel(bc.face);
    if numel(opt.(name)) == 1
        opt.(name) = repmat(opt.(name), nf, 1);
    end
    assert(numel(opt.(name)) == nf, ...
                  [name, ' must be a scalar, or one per face in bc.face']);
end