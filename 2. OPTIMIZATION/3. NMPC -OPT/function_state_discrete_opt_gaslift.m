function xk1 = function_state_discrete_opt_gaslift(xk, uk)
% Input:
% Qginj = u(1)                Injected gas flow rate (cm3/s)
                        
% States: 
% mg = x(1)                   Gass mass in the pipe (g)
% ml = x(2)                   Liquid mass in the pipe (g)
%      x(3)                   State of unmeasured disturbance

% Output: 
% wlout = y(1)                Liquid outflow rate (g/s)

%% Specify parameters
sample_time = 0.05;
N = 10;                       % Number of integration time steps for Euler method
dt = sample_time / N;         % Step size

%% Multistep forward Euler integration method
uk1 = [uk; 0];
xk1 = xk(:);

for i = 1:N
    xk1 = xk1 + dt * function_state_continuous_opt_gaslift(xk1, uk1);
end

