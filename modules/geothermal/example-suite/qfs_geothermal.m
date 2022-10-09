function [description, options, state0, model, schedule, plotOptions] = qfs_geothermal(varargin)

    options = struct('NaCl', false);
    [options, extra] = merge_options(options, varargin{:});
    [~, options_qfs_wo, state0, model, schedule, plotOptions] = qfs_wo(extra{:});
    
    options = merge_structs(options, options_qfs_wo);
    description = '';
    if nargout <= 2, return; end
    
    G     = model.G;
    rock  = model.rock;
    fluid = initSimpleADIFluid('phases', 'W'                  , ...
                               'n'     , 1                    , ...
                               'rho'   , 1000*kilogram/meter^3, ...
                               'mu'    , 0.1*centi*poise      , ...
                               'c'     , 1e-10/Pascal         , ...
                               'pRef'  , 1*atm                );
    
    fluid = addThermalFluidProps(fluid, 'useEOS', false               , ...
                                        'cT'    , 1e-6/Kelvin         , ...
                                        'cX'    , log(2)              , ...
                                        'TRef'  , (273.15 + 20)*Kelvin);
    
    rock  = addThermalRockProps(rock);
    
    compFluid = [];
    if options.NaCl
        compFluid = CompositionalBrineFluid(          ...
            {'H2O'             , 'NaCl'            }, ... % Names
            [18.015281*gram/mol, 58.442800*gram/mol], ... % Molar masses
            [0                 , 1e-6              ]);    % Molecular diffusivities
    end
    model = GeothermalModel(G, rock, fluid, compFluid);
    
    state0.s = ones(G.cells.num, 1);
    state0.T = repmat((273.15 + 20)*Kelvin, G.cells.num, 1);

    state0.components = ones(G.cells.num,1);
    if options.NaCl
        x0 = repmat(0.1, G.cells.num, 1);
        state0.components = [state0.components - x0, x0];
    end
    
    [schedule.control(1).W] = deal(addThermalWellProps(schedule.control(1).W, 'T', 273.15*Kelvin + 200));
    [schedule.control(1).W.components] = deal(1);
    [schedule.control.W.compi] = deal(1);
    
end