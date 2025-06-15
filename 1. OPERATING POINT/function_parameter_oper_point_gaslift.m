function  par = function_parameter_oper_point_gaslift
%% Parameters
par.D = 2.54;                                      % Diameter of the vertical pipe (cm)
par.A = pi * par.D^2 / 4;                          % Area of the vertical pipe (cm2)
par.g = 980.665;                                   % Acceleration of gravity (cm/s2)
par.H = 426.72;                                    % Height of the vertical pipe (cm)
par.S = 0.532;                                     % Ratio of the height of the vertical pipe to height of the reservoir
par.Hres = par.H * par.S;                          % Height of the reservoir (cm)
par.V = par.A * par.H;                             % Volume of the vertical pipe (cm3)
par.Mg = 28.951;                                   % Molecular weight of gas (g/mol)
par.patm = 1013250;                                % Atmospheric pressure (g/cm-s2 (1e-1*Pa))
par.R = 8.314462618 * 1e7;                         % Universal gas constant (cm3⋅(g/cm-s2)/(K⋅mol))
par.rhol = 0.998204;                               % Liquid density (g/cm3) at 20 °C (Perry's Chemical Engineering Handbook, page 2-91/139)
par.pres = par.patm + par.rhol * par.g * par.Hres; % Reservoir pressure (hydrostatic pressure) (63.5 cm H2O = 6227.223 Pa) (g/cm-s2 = 1e-1*Pa)
par.T = 293.15;                                    % Temperature at the vertical pipe top (K) (20 °C)

%% Initial conditions
par.alphal0 = 0.532;                               % Initial liquid volume fraction
par.ml0 = par.V * par.alphal0 * par.rhol;          % Initial liquid mass (g)
par.ptop0 = 1.01 * par.patm;                       % Initial pressure at the top of the pipe (g/cm-s2 (1e-1*Pa))
par.mg0 = par.ptop0 * par.V * (1 - par.alphal0) * par.Mg / (par.R * par.T); % Initial gas mass (g)

function_align_comments('function_parameter_oper_point_gaslift.m', 51)
