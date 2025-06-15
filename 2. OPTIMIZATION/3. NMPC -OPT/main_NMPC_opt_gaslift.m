clear all
clc
tic
%% Linear time-invariant state-space MPC of a gas-lift system at equilibrium 
% Description: 
% A vertical pipe is submerged in a water-filled reservoir. When the air is 
% injected into the bottom of the pipe, water starts to flow out. However,  
% without air injection, there is no water outflow from the pipe. 
%
% Input:
% Qginj = u(1)    Injected gas flow rate (cm3/s) [manipulated variable]
%
% States: 
% mg = x(1)       Gass mass in the pipe (g)
% ml = x(2)       Liquid mass in the pipe (g)
%      x(3)       State of unmeasured disturbance
%
% Output: 
% wlout = y(1)    Liquid outflow rate (g/s)

%% Get parameters.
par = function_parameters_opt_gaslift;

% % Store parameters in an Excel file.
% T = {'Run', 'Sample Time (s)', 'Prediction Horizon Steps', 'Control Horizon Steps', 'Comment'};
% writecell(T, 'parameters_NMPC_opt_gaslift.xlsx', 'Sheet', 1, 'Range', 'A1:E1');
% T = table({'1'}, par.sample_time, par.prediction_horizon, par.control_horizon);
% writetable(T, 'parameters_NMPC_opt_gaslift.xlsx', 'Sheet', 1, 'WriteVariableNames', 0, 'Range', 'A2:D2');
% winopen('parameters_NMPC_opt_gaslift.xlsx')

%% Configure nonlinear MPC.
% nlobj = nlmpc(nx,ny,'MV',mvIndex,'MD',mdIndex,'UD',udIndex)
nlobj = nlmpc(3, 1,'MV', 1, 'UD', 2); 

% Specify horizon:
% The prediction model sample time is the same as the controller sample time.
nlobj.Ts = par.sample_time;
nlobj.PredictionHorizon = par.prediction_horizon; 
nlobj.ControlHorizon = par.control_horizon;

% Specify state and output functions.
nlobj.Model.StateFcn = 'function_state_continuous_opt_gaslift';
nlobj.Model.OutputFcn = 'function_output_opt_gaslift';

% Specify cost function. 
nlobj.Optimization.CustomCostFcn = "function_cost_opt_end_gaslift";
nlobj.Optimization.ReplaceStandardCost = true;

% % Set minimum and maximum manipulated variable rates (injected gas flow rate [m3/s]).
% nlobj.ManipulatedVariables.RateMin = -100;
% nlobj.ManipulatedVariables.RateMax = 100;
% 
% % Set a weight for the manipulated variable rate.
% nlobj.Weights.ManipulatedVariablesRate = 1;

% Specify a lower bound for the manipulated variable.
nlobj.ManipulatedVariables.Min = 0;

%  Specify initial conditions.
x0 = [par.mg0; par.ml0; 0];   % States (g)
u0 = par.Qginj0;              % Inputs (cm3/s)
%validateFcns(nlobj, x0(1:3), par.Qginj0);

%% Open and simulate the model.
model = 'NMPC_opt_gaslift';
open_system(model)
sim(model)

%% Plot figures.
fig = figure;

subplot(1, 2, 1)
plot(logsout{18}.Values.Time, logsout{18}.Values.Data, '-r', 'LineWidth', 1) % Injected gas flow rate (cm3/s)
xlabel('Time (s)')
ylabel('Q_{g, inj} (cm^3/s) (MV)', 'Interpreter', 'tex')
title('')

subplot(1, 2, 2)
plot(logsout{17}.Values.Time, logsout{17}.Values.Data, '-r', 'LineWidth', 1) % Liquid outflow rate (g/s)
xlabel('Time (s)')
ylabel('w_{l, out} (g/s) (OV)', 'Interpreter', 'tex')
legend('ph = 200, ch = 40, T_s = 0.05',  'Location', 'northoutside', 'Orientation','horizontal')
legend('boxoff')
title('')
% savefig(fig, 'fig5_NMPC_opt_ph200_ch40_Ts0_05_Q1e-3_gaslift.fig') 
% close(fig)
toc
