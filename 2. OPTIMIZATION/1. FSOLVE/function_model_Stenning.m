function [F, wlout] = function_model_Stenning(x, Qginj, par)
%% Provides the liquid outflow rate for the given gas inflow rate at steady-state
% Input: air inflow rate, Qginj (cm3/s)
% Output: water outflow rate, wlout (g/s)

%% States
mg = x(1); % Gas mass (g)
ml = x(2); % Liquid mass (g)

%% Algebraic equations for the vertical pipe
% Density (g/cm3)
    % Volume of liquid (cm3)
    Vl = ml / par.rhol;
    % Average gas density 
    rhog_avg = mg / (par.V - Vl);
    % Average mixture density
    rhomix_avg = (mg + ml) / par.V;
% Pressure (g/cm-s2)
    % Hydrostatic pressure
    ph = rhomix_avg * par.g * par.H;
    % Average pressure
    pavg = rhog_avg * par.R * par.T / par.Mg;
    % Pressure at the top of the pipe
    ptop = pavg - ph / 2;
    % Gas density at the top of the pipe (g/cm3)
    rhogtop = ptop * par.Mg / (par.R * par.T);
    % Pressure at the pipe bottom 
    pbot = ptop + ph;
% Volume fraction
    % Average liquid volume fraction 
    alphal_avg = Vl / par.V;
    % Average gas volume fraction
    alphag_avg = 1 - alphal_avg;
    % Gas volume fraction at the pipe top
    alphagtop = alphag_avg * pavg/ptop;
    % Gas volume fraction at the pipe bottom 
    alphagbot = 2 * alphag_avg - alphagtop;
% Density (g/cm3)
    % Gas injection density, injection pressure is average pressure in the pump
    rhoginj = pbot * par.Mg / (par.R * par.T); 
% Liquid flow rate (g/s)
    % Liquid inflow rate 
    wlin = (1 - alphagbot) * par.A * sqrt(par.rhol * max(0, par.pres - pbot));
    % Liquid outflow rate 
    wlout = (1 - alphagtop) * par.A * sqrt(par.rhol * max(0, ptop - par.patm));
% Volumetric flow rate (cm3/s)
    % Liquid volumetric outflow rate 
    Qlout = wlout / par.rhol;
    % Gas volumetric outflow rate
    Qgout = Qlout * (alphagtop) / ((1 - alphagtop)^2);
% Gas flow rate 
    % Gas mass inflow rate (g/s)
    wginj = Qginj * rhoginj;
    % Gas mass outflow rate (g/s)
    wgout = Qgout * rhogtop;

%% State equations
F(1) = wginj - wgout; % Rate of change of gas mass (g/s)
F(2) = wlin - wlout; % Rate of change of liquid mass (g/s)