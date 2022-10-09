function [ mu_b ] = viscosity_brine( p,t,c )
% Compute viscosity of a brine following Spivey et al., 2004
%
% SYNOPSIS:
%   mu_b = viscosity_brine(p,t,c)
% 
% PARAMETERS: 
%   p   - pressure. One value per cell. 
%   t   - temperature. One value per cell. 
%   c   - NaCl mass fraction. One value per cell.
% 
% RETURNS:
%   mu_b - brine viscosity. One value per cell 
% 
% SEE ALSO:
%   'density_pure_water', 'viscosity_pure_water', 'density_brine'
p               = p./1e6;                                                   % conversion from Pa tp MPa 
t               = t-273.15;                                                 % conversion from Kelvin to �Celsius
MW_NaCl         = 0.0584428;                                                % [kg/mol] Molar weigth of NaCl
X_NaCl          = c;                                                        % mass fraction of NaCl  
X_Water         = 1 - X_NaCl;                                               % mass fraction of water    
cNaCl           = X_NaCl./X_Water;                                          % concentration of NaCl per kg of water [kg/kg]   
mol_NaCl        = cNaCl./MW_NaCl;                                           % molality [mol/kg H2O]  

tMin            = 1e-3; % Greater than 0 to avoid singular term1, term2
tMax            = 275;
pMin            = 0;
pMax            = 200;
mol_NaClMin     = 0;
mol_NaClMax     = 5.7;

if max(value(t)) > 275
    warning('T out of range')
end

if max(value(p)) > 200
    warning('p out of range')
end

if max(value(mol_NaClMax)) > 5.7
    warning('mol_NaCl out of range')
end

%% Caping for the ADI variable scheme
t = min(max(t, tMin), tMax);
p = min(max(p, pMin), pMax);
mol_NaCl = min(max(mol_NaCl, mol_NaClMin), mol_NaClMax);


%% Coefficient for the scaling
A1 = 0.0173; 
A2 = 0.068;
B1 = -1.0531; 
B2 = 0.0273;


mu_kestin = kestin_brine_viscosity(p, 125,mol_NaCl); % function kestin_brine_viscosity take MPa in input
mu_k      = kestin_brine_viscosity(p,t,mol_NaCl);

term1 = (A2.*(p./100)+A1).*(log(t./125)).^2;
term2 = (B2.*(p./100)+B1).*(log(t./125));

mu    = mu_kestin.*exp(term1 + term2);

% if value(t) > 125 
%     mu_b = mu;
% else 
%     mu_b = mu_k;
% end
%     mu_b = mu_b.*1e-6; 

useKestin = value(t) < 125;
mu_b = mu.*(~useKestin) + mu_k.*(useKestin);
mu_b = mu_b.*1e-6;

end

%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

%%--------------------------------------------------------------------------
%%                             References
%%--------------------------------------------------------------------------
% Spivey, J.P., Mccain, W.D., North, R. 2004. Estimating density, formation volume factor, compressibility,
% methane solubility, and viscosity for oilfield brines at temperatures from 0 to 275ºC, pressures to 200 Mpa,
% and salinities to 5.7 mol.kg-1, Journal of Canadian Petroleum Technology, 43 (7), 52-60.
