function  par = function_parameters_NMPC_PMM_mg_gaslift
%% Control parameters
% Reference for the liquid outflow rate (g/s)
par.wloutref = 300; 

% Horizon steps
par.prediction_horizon = 10;             
par.control_horizon = 2;                

% Simulation parameters (s)
par.sample_time = 0.1; 
par.simulation_time = 100; 

%% System parameters
% Input parameters
par.D = 2.54;                 % Diameter of the vertical pipe (cm)
par.H = 426.72;               % Height of the vertical pipe (cm)
par.S = 0.532;                % Ratio of the height of the vertical pipe to height of the reservoir
par.rhol = 0.998204;          % Liquid density (g/cm3) at 20 °C
par.T = 293.15;               % Temperature (K) [20 °C]

% Constants
par.g = 980.665;              % Acceleration of gravity (cm/s2)
par.Mg = 28.951;              % Molecular weight of gas (g/mol)
par.patm = 1013250;           % Atmospheric pressure (g/cm-s2 [0.1 * Pa])
par.R = 8.314462618 * 1e7;    % Universal gas constant (g-cm2/s2-K-mol)

% Parameter calculations
par.A = pi * par.D^2 / 4;     % Area of the vertical pipe (cm2)
par.Hres = par.H * par.S;     % Height of the reservoir (cm)
par.pres = par.patm + par.rhol * par.g * par.Hres; % Reservoir pressure (hydrostatic pressure) (g/cm-s2 [0.1 * Pa])
par.V = par.A * par.H;        % Volume of the vertical pipe (cm3)

%% References for physical and chemical parameters
% <https://physics.nist.gov/cuu/Constants/index.html> 
% <https://webbook.nist.gov/chemistry/fluid/> 

%% Initial conditions
% Initial input
par.Qginj0 = 9.43152759202956; % Injected gas flow rate [MV] (cm3/s)

% Initial states
par.mg0 = 0.135328793583754;  % Gas mass (g)
par.ml0 = 1146.84497366167;   % Liquid mass (g)
 
% function_align_comments('function_parameters_NMPC_PMM_mg_gaslift.m', 30)
