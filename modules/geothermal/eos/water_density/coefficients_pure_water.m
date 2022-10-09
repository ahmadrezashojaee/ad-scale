%% Fonction that computes the coefficients for equation of state of pure water
% returns the coefficients Ew and Fw as well as the density of pure water
% as a function of temperature for the reference pressure p0 = 70 Mpa

function [ rho_w0, Ew, Fw ] = coefficients_pure_water( t )

T           = t/100;

%% Table 2, Spivey et al., 2004
Dw_v        = [-0.127213, 0.645486, 1.03265, -0.070291, 0.639589];
Ew_v        = [4.221, -3.478, 6.221, 0.5182, -0.4405];
Fw_v        = [-11.403, 29.932, 27.952, 0.20684, 0.3768];


%% Compute coefficients from Eq.(3) and (9) in Spivey et al., 2004
rho_w0 = (Dw_v(1).*(T.^2) + Dw_v(2).*T + Dw_v(3)) ./ (Dw_v(4).*(T.^2) + Dw_v(5).*T +1);
Ew     = (Ew_v(1).*(T.^2) + Ew_v(2).*T + Ew_v(3)) ./ (Ew_v(4).*(T.^2) + Ew_v(5).*T +1);
Fw     = (Fw_v(1).*(T.^2) + Fw_v(2).*T + Fw_v(3)) ./ (Fw_v(4).*(T.^2) + Fw_v(5).*T +1);

end

