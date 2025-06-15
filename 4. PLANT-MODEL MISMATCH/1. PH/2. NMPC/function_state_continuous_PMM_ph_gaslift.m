function dxdt = function_state_continuous_PMM_ph_gaslift(x, u) 
%% Computes the rate of change of the states.

% Output 
    % wlout = y(1)            Liquid outflow rate (g/s)

%% Input
Qginj = u(1);                 % Injected gas flow rate (cm3/s)
%       u(2)                  Rate of change of unmeasured disturbance

%% States 
mg = x(1);                    % Gass mass in the pipe (g)
ml = x(2);                    % Liquid mass in the pipe (g)
%    x(3)                     % Unmeasured disturbance

%% System parameters
% Input parameters
D = 2.54;                     % Diameter of the vertical pipe (cm)
H = 426.72;                   % Height of the vertical pipe (cm)
S = 0.532;                    % Ratio of the height of the vertical pipe to height of the reservoir
rhol = 0.998204;              % Liquid density (g/cm3) at 20 °C 
T = 293.15;                   % Temperature (K) [20 °C]

% Constants
g = 980.665;                  % Acceleration of gravity (cm/s2)
Mg = 28.951;                  % Molecular weight of gas (g/mol)
patm = 1013250;               % Atmospheric pressure (g/cm-s2 [0.1 * Pa])
R = 8.314462618 * 1e7;        % Universal gas constant (g-cm2/s2-K-mol)

% Parameter calculations
A = pi * D^2 / 4;             % Area of the vertical pipe (cm2)
Hres = H * S;                 % Height of the reservoir (cm)
pres = patm + rhol * g * Hres;% Reservoir pressure (hydrostatic pressure) (g/cm-s2 [0.1 * Pa])
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
    % Hydrostatic pressure --> ph = rhomix_avg * g * H changed to
    ph = rhomix_avg * H;
    % Average pressure
    pavg = rhog_avg * R * T / Mg;
    % Pressure at the top of the pipe
    ptop = pavg - ph / 2;
    % Pressure at the pipe bottom 
    pbot = ptop + ph;
% Gas density at the top of the pipe (g/cm3)
    rhogtop = ptop * Mg / (R * T);
% Volume fraction
    % Average liquid volume fraction 
    alphal_avg = Vl / V;
    % Average gas volume fraction
    alphag_avg = 1 - alphal_avg;
    % Gas volume fraction at the pipe top
    alphagtop = alphag_avg * pavg / ptop;
    % Gas volume fraction at the pipe bottom 
    alphagbot = 2 * alphag_avg - alphagtop;
% Liquid flow rate (g/s)
    % Liquid inflow rate --> Signed square root
    a = rhol * (pres - pbot);
    wlin = (1 - alphagbot) * A * sign(a) * sqrt(abs(a));   
    % Liquid outflow rate --> Signed square root
    b = rhol * (ptop - patm);
    wlout = (1 - alphagtop) * A * sign(b)* sqrt(abs(b));
% Volumetric flow rate (cm3/s) 
    % Liquid volumetric outflow rate 
    Qlout = wlout / rhol;
    % Gas volumetric outflow rate 
    Qgout = Qlout * (alphagtop) / ((1 - alphagtop)^2);
% Gas injection density (g/cm3) --- injection pressure is the bottom hole pressure
    rhoginj = pbot * Mg / (R * T); 
% Gas flow rate (g/s)
    % Gas mass inflow rate 
    wginj = Qginj * rhoginj;
    % Gas mass outflow rate 
    wgout = Qgout * rhogtop;

%% State equations
dxdt = zeros(3, 1);
dxdt(1) = wginj - wgout;      % Rate of change of gas (g/s)
dxdt(2) = wlin - wlout;       % Rate of change of liquid (g/s)
dxdt(3) = u(2);               % Rate of change of unmeasured disturbance



