clear all
clc
% tic
%% Linear time-invariant state-space MPC of a gas-lift system at equilibrium 
% Description: 
% A vertical pipe is submerged in a water-filled reservoir. When the air is 
% injected into the bottom of the pipe, water starts to flow out. However,  
% without air injection, there is no water outflow from the pipe. 
%
% Input:
% Qginj = u(1)    Injected gas flow rate (cm3/s) [manipulated variable]
% pres  = u(2)    Reservoir pressure (g/cm-s2) [measured disturbance]
%
% States: 
% mg = x(1)       Gass mass in the pipe (g)
% ml = x(2)       Liquid mass in the pipe (g)
%      x(3)       State of unmeasured disturbance
%
% Output: 
% wlout = y(1)    Liquid outflow rate (g/s)

%% Get parameters.
par = function_parameters_NMPC_MD_gaslift;

% % Store parameters in an Excel file.
% T = {'Run', 'Sample Time (s)', 'Reservoir Pressure Change [g/s])', 'Comment'};
% writecell(T, 'parameters_NMPC_MD_gaslift.xlsx', 'Sheet', 1, 'Range', 'A1:D1');
% T = table({'1'}, par.sample_time, par.wloutref);
% writetable(T, 'parameters_NMPC_MD_gaslift.xlsx', 'Sheet', 1, 'WriteVariableNames', 0, 'Range', 'A2:C2');
% winopen('parameters_NMPC_MD_gaslift.xlsx')

%% Configure nonlinear MPC. 
% nlobj = nlmpc(nx,ny,'MV',mvIndex,'MD',mdIndex,'UD',udIndex)
nlobj = nlmpc(3, 1,'MV', 1, 'MD', 2, 'UD', 3); 

% Specify horizon:
% The prediction model sample time is the same as the controller sample time.
nlobj.Ts = par.sample_time;
nlobj.PredictionHorizon = par.prediction_horizon; 
nlobj.ControlHorizon = par.control_horizon;

% Specify state and output functions.
nlobj.Model.StateFcn = 'function_state_continuous_MD_gaslift';
nlobj.Model.OutputFcn = 'function_output_MD_gaslift';

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
u0 = [par.Qginj0; par.pres];  % Inputs (cm3/s, g/cm-s2 [0.1 * Pa])
%validateFcns(nlobj, x0(1:3), par.Qginj0, par.pres);

%% Open and simulate the model.
model = 'NMPC_MD_gaslift';
open_system(model)
sim(model)

%% Plot figures.
fig = figure;

subplot(1, 2, 1)
plot(logsout{19}.Values.Time, logsout{19}.Values.Data, '-r', 'LineWidth', 1) % Injected gas flow rate (cm3/s)
xlabel('Time (s)')
ylabel('Q_{g, inj} (cm^3/s) (MV)', 'Interpreter', 'tex')
legend('NMPC', 'Location', 'northoutside')
legend('boxoff')
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

% savefig(fig, 'fig2_NMPC_MD_0_02pres_Ts0_1_ph10_ch2_gaslift.fig') 
% close(fig)

