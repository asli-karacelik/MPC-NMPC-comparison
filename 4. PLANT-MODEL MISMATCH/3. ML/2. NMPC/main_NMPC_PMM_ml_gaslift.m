clear all
clc
%% Linear time-invariant state-space model predictive control (MPC) for a gas-lift system at equilibrium
% Description: 
    % A vertical pipe is submerged in a water-filled reservoir. When the air is 
    % injected into the bottom of the pipe, water starts to flow out. However,  
    % without air injection, there is no water outflow from the pipe. 

% Input
    % Qginj = u(1)            Injected gas flow rate (cm3/s) [manipulated variable]
    %         u(2)            Rate of change of unmeasured disturbance
    
% States 
    % mg = x(1)               Gass mass in the pipe (g)
    % ml = x(2)               Liquid mass in the pipe (g)
    %      x(3)               Unmeasured disturbance

% Output 
    % wlout = y(1)            Liquid outflow rate (g/s)

%% Get parameters.
par = function_parameters_NMPC_PMM_ml_gaslift;

% % Store parameters in an Excel file
% % T = {'Run', 'Sample Time (s)', 'Prediction Horizon', 'Control Horizon', 'Comment'};
% % writecell(T, 'parameters_NMPC_PMM_ml_gaslift.xlsx', 'Sheet', 1, 'Range', 'A1:D1');
% T = table({'4'}, par.sample_time, par.prediction_horizon, par.control_horizon);
% writetable(T, 'parameters_NMPC_PMM_ml_gaslift.xlsx', 'Sheet', 1, 'WriteVariableNames', 0, 'Range', 'A5:D5');
% winopen('parameters_NMPC_PMM_ml_gaslift.xlsx')

%% Configure nonlinear NMPC. 
% nlobj = nlmpc(nx,ny,'MV',mvIndex,'MD',mdIndex,'UD',udIndex)
nlobj = nlmpc(3, 1,'MV', 1, 'UD', 2); 

% The prediction model sample time is the same as the controller sample time.
nlobj.Ts = par.sample_time;

% Specify horizon.
nlobj.PredictionHorizon = par.prediction_horizon; 
nlobj.ControlHorizon = par.control_horizon;

% Specify state and output functions.
nlobj.Model.StateFcn = 'function_state_continuous_PMM_ml_gaslift';
nlobj.Model.OutputFcn = 'function_output_PMM_ml_gaslift';

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
% validateFcns(nlobj, x0(1:3), par.Qginj0);

%% Open and simulate the model.
model = 'NMPC_PMM_ml_gaslift';
open_system(model)
sim(model)

%% Plot figures.
fig = figure;

subplot(1, 2, 1)
plot(logsout{19}.Values.Time, logsout{19}.Values.Data, '-r', 'LineWidth', 1) % Injected gas flow rate (cm3/s)
xlabel('Time (s)')
ylabel('Q_{g, inj} (cm^3/s) (MV)', 'Interpreter', 'tex')
title('')

subplot(1, 2, 2)
plot(logsout{1}.Values.Time, logsout{1}.Values.Data, '-k', 'LineWidth', 1) % Reference for liquid outflow rate (g/s)
hold on
plot(logsout{18}.Values.Time, logsout{18}.Values.Data, '-r', 'LineWidth', 1) % Liquid outflow rate (g/s)
xlabel('Time (s)')
ylabel('w_{l, out} (g/s) (OV)', 'Interpreter', 'tex')
legend('Reference', 'NMPC',  'Location', 'northoutside', 'Orientation','horizontal')
legend('boxoff')
title('')

% savefig(fig, 'fig5_MPC_PMM_ml_Ts0_05_ph10_ch2_gaslift.fig') 
% close(fig)

