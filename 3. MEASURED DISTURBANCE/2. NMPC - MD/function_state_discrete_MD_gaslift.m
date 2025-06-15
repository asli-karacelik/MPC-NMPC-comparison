function xk1 = function_state_discrete_MD_gaslift(xk, uk) 
% Input:
% Qginj = u(1)                Injected gas flow rate (cm3/s) [manipulated variable]
% pres = u(2)                 Reservoir pressure (g/cm-s2) [measured disturbance]
%
% States: 
% mg = x(1)                   Gass mass in the pipe (g)
% ml = x(2)                   Liquid mass in the pipe (g)
%      x(3)                   State of unmeasured disturbance
%
% Output: 
% wlout = y(1)                Liquid outflow rate (g/s)

%% Specify parameters
sample_time = 0.05;
N = 10;                       % Number of integration time steps for Euler method
dt = sample_time / N;         % Step size

%% Multistep forward Euler integration method
uk1 = [uk(:); 0];
xk1 = xk(:);

for i = 1:N
    xk1 = xk1 + dt * function_state_continuous_MD_gaslift(xk1, uk1);
end

