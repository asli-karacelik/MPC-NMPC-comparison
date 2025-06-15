function y = function_measurement_opt_gaslift(x)  % Output is the liquid outflow rate (kg/s)
% Input:
% Qginj = u(1)                % Injected gas flow rate (cm3/s)

% States: 
mg = x(1);                    % Gass mass in the pipe (g)
ml = x(2);                    % Liquid mass in the pipe (g)
%    x(3)                     % State of unmeasured disturbance

% Output: 
% wlout = y(1)                % Liquid outflow rate (g/s)

%% System parameters
% Input parameters
D = 2.54;                     % Diameter of the vertical pipe (cm)
H = 426.72;                   % Height of the vertical pipe (cm)
rhol = 0.998204;              % Liquid density (g/cm3) at 20 °C 
T = 293.15;                   % Temperature (K) [20 °C]

% Constants
g = 980.665;                  % Acceleration of gravity (cm/s2)
Mg = 28.951;                  % Molecular weight of gas (g/mol)
patm = 1013250;               % Atmospheric pressure (g/cm-s2 [0.1 * Pa])
R = 8.314462618 * 1e7;        % Universal gas constant (g-cm2/s2-K-mol)

% Parameter calculations
A = pi * D^2 / 4;             % Area of the vertical pipe (cm2)
V = A * H;                    % Volume of the vertical pipe (cm3)

%% References for physical and chemical parameters
% <https://physics.nist.gov/cuu/Constants/index.html> 
% <https://webbook.nist.gov/chemistry/fluid/> 

%% Algebraic equations 
% Volume of liquid (cm3)
    Vl = ml / rhol;
% Density (g/cm3)
    % Average gas density 
    rhog_avg = mg / (V - Vl);
    % Average mixture density
    rhomix_avg = (mg + ml) / V;
% Pressure (g/cm-s2)
    % Hydrostatic pressure
    ph = rhomix_avg * g * H;
    % Average pressure
    pavg = rhog_avg * R * T / Mg;
    % Pressure at the top of the pipe
    ptop = pavg - ph / 2;
% Volume fraction
    % Average liquid volume fraction 
    alphal_avg = Vl / V;
    % Average gas volume fraction
    alphag_avg = 1 - alphal_avg;
    % Gas volume fraction at the pipe top
    alphagtop = alphag_avg * pavg/ptop;
% Liquid outflow rate (g/s)
    wlout = (1 - alphagtop) * A * sqrt(rhol * max(0, ptop - patm));

%% Measurement function
y = wlout + x(3);
