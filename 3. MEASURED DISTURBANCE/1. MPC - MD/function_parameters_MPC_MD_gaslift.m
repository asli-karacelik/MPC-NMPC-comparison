function  par = function_parameters_MPC_MD_gaslift
%% Control parameters
% Reference for the liquid outflow rate (g/s)
par.wloutref = 300; 

% Horizon steps
par.prediction_horizon = 50;             
par.control_horizon = 10;                

% Simulation parameters (s)
par.sample_time = 0.1; 
par.simulation_time = 300; 

%% System parameters
% Input parameters
par.D = 2.54;                             % Diameter of the vertical pipe (cm)
par.H = 426.72;                           % Height of the vertical pipe (cm)
par.S = 0.532;                            % Ratio of the height of the vertical pipe to height of the reservoir
par.rhol = 0.998204;                      % Liquid density (g/cm3) at 20 °C 
par.T = 293.15;                           % Temperature (K) [20 °C]

% Constants
par.g = 980.665;                          % Acceleration of gravity (cm/s2)
par.Mg = 28.951;                          % Molecular weight of gas (g/mol)
par.patm = 1013250;                       % Atmospheric pressure (g/cm-s2 [0.1 * Pa])
par.R = 8.314462618 * 1e7;                % Universal gas constant (g-cm2/s2-K-mol)

% Parameter calculations
par.A = pi * par.D^2 / 4;                 % Area of the vertical pipe (cm2)
par.Hres = par.H * par.S;                 % Height of the reservoir (cm)
par.V = par.A * par.H;                    % Volume of the vertical pipe (cm3)

%% Measured disturbance: reservoir pressure (hydrostatic pressure) (g/cm-s2 [0.1 * Pa])
par.pres = par.patm + par.rhol * par.g * par.Hres;   
par.pres_final = 0.98 * par.pres;         

%% References for physical and chemical parameters
% <https://physics.nist.gov/cuu/Constants/index.html> 
% <https://webbook.nist.gov/chemistry/fluid/> 

%% Initial conditions
par.alphal0 = 0.532;                 % Initial liquid volume fraction
par.ptop0 = 1.01 * par.patm;         % Initial top pressure (g/cm-s2 [0.1 * Pa])

% Initial states
par.ml0 = par.V * par.alphal0 * par.rhol; % Initial liquid mass (g)
par.mg0 = par.ptop0 * par.V * (1 - par.alphal0) * par.Mg / (par.R * par.T); % Initial gas mass (g)

%function_align_comments('function_parameters_MPC_MD_gaslift.m', 40)
