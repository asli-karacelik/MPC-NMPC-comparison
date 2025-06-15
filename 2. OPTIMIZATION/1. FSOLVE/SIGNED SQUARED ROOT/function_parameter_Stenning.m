function par = function_parameter_Stenning(S)
%% Parameters
% Input parameters
par.S = S;                           % Submergence ratio: ratio of the height of the vertical pipe to height of the reservoir
par.D = 2.54;                        % Diameter of the vertical pipe (cm)
par.H = 426.72;                      % Height of the vertical pipe (cm)
par.rhol = 0.99821;                  % Liquid density (g/cm3) at 20 °C
par.T = 293.15;                      % Temperature (K) [20 °C]

% Constants
par.g = 980.665;                     % Acceleration of gravity (cm/s2)
par.Mg = 28.9586;                    % Molecular weight of air (g/mol)
par.patm = 1013250;                  % Atmospheric pressure (g/cm-s2 [0.1 * Pa])
par.R = 8.314462618 * 1e7;           % Universal gas constant (g-cm2/s2-K-mol)

% Parameter calculations
par.A = pi * par.D^2 / 4;            % Area of the vertical pipe (cm2)
par.Hres = par.H * par.S;            % Height of the reservoir (cm)
par.V = par.A * par.H;               % Volume of the vertical pipe (cm3)
par.pres = par.patm + par.rhol * par.g * par.Hres; % Reservoir pressure (g/cm-s2 [0.1 * Pa])

%% Initial conditions
par.alphal0 = 0.1;                   % Initial liquid volume fraction
par.ptop0 = par.patm;                % Initial top pressure (g/cm-s2 [0.1 * Pa])

% Initial states
par.ml0 = par.V * par.alphal0 * par.rhol; % Initial liquid mass (g)
par.mg0 = par.ptop0 * par.V * (1 - par.alphal0) * par.Mg / (par.R * par.T); % Initial gas mass (g)

end

% align_comments('function_parameter_Stenning.m', 40)