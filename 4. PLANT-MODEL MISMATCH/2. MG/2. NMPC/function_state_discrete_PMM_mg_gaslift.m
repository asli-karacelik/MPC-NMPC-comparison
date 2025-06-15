function xk1 = function_state_discrete_PMM_mg_gaslift(xk, uk) 
%% Implements a multi-step forward Euler integration (explicit) method for updating the states of a gas-lift system.

% Input
    % Qginj = uk(1)           Injected gas flow rate (cm3/s) [manipulated variable]
    %         uk(2)           Rate of change of unmeasured disturbance
    
% States 
    % mg = xk(1)              Gass mass in the pipe (g)
    % ml = xk(2)              Liquid mass in the pipe (g)
    %      xk(3)              Unmeasured disturbance

% Outputs
    % xk1                     Next states after discrete-time update 

%% Specify parameters
sample_time = 0.1;
N = 100;                       % Number of integration time steps for Euler method
dt = sample_time / N;         % Step size

%% Initilization
% Specify current inputs and states.
uk = [uk; 0];
xk = xk(:); 
% Initialize next states
xk1 = xk;

%% Multi-step forward Euler integration
for i = 1:N
    % Compute state derivatives at the current states and inputs.
    dxdt = function_state_continuous_PMM_mg_gaslift(xk1, uk);
    
    % Update states using Euler integration
    xk1 = xk1 + dt * dxdt;
end


