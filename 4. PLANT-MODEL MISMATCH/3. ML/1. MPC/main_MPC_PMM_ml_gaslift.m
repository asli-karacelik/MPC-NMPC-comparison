clear all
clc
%% Linear time-invariant state-space MPC of a gas-lift system at equilibrium 
% Description: 
%   A vertical pipe is submerged in a water-filled reservoir. When the air is 
%   injected into the bottom of the pipe, water starts to flow out. However,  
%   without air injection, there is no water outflow from the pipe. 
%
% Input:
%   - Qginj = u(1)    Injected gas flow rate (cm3/s) [manipulated variable]
%
% States: 
%   - mg = x(1)       Gass mass in the pipe (g)
%   - ml = x(2)       Liquid mass in the pipe (g)
%   -      x(3)       Unmeasured disturbance
%
% Output: 
%   - wlout = y(1)    Liquid outflow rate (g/s)

%% Get parameters.
par = function_parameters_MPC_PMM_ml_gaslift;

% % Store parameters in an Excel file
% T = {'Run', 'Sample Time (s)', 'Prediction Horizon', 'Control Horizon', 'Covariance Matrix (Q)', 'Comment'};
% writecell(T, 'parameters_MPC_PMM_ml_gaslift.xlsx', 'Sheet', 1, 'Range', 'A1:F1');
% T = table({'1'}, par.sample_time, par.prediction_horizon, par.control_horizon);
% writetable(T, 'parameters_MPC_PMM_ml_gaslift.xlsx', 'Sheet', 1, 'WriteVariableNames', 0, 'Range', 'A2:D2');
% % winopen('parameters_MPC_PMM_gaslift.xlsx')

%% Linearize the plant model. 
% The internal plant model is obtained by linearizing a nonlinear plant model,
% |model_PMM_gaslift|, in Simulink at a nonzero steady-state operating point.

% Generate an operating point specification for the model's initial conditions.
plant_model = 'model_PMM_ml_gaslift';
sim(plant_model)
operating_point_specification = operspec(plant_model);

% Obtain the operating point from the model model's initial conditions.
[operating_point, operating_report] = findop(plant_model, operating_point_specification);

% Specify nominal state, output, and input values from the operating point.
x0 = [operating_report.States(1).x; operating_report.States(2).x];
y0 = operating_report.Outputs.y;
u0 = operating_report.Inputs.u;

% Linearize the plant at the operating point.
plant = linearize(plant_model, operating_point); 

%% Configure MPC. 
% mpcobj = mpc(plant, ts, P, M, W, MV, OV, DV)
mpcobj = mpc(plant, par.sample_time, par.prediction_horizon, par.control_horizon);

% Set the nominal values in the controller.
mpcobj.Model.Nominal = struct('X', x0, 'U', u0, 'Y', y0);

% Change the plant model from continuous to discrete time.
plant_discrete = c2d(plant, par.sample_time);

% % Set minimum and maximum manipulated variable rates (injected gas flow rate [m3/s]).
% mpcobj.ManipulatedVariables = struct('RateMin', -100, 'RateMax', 100); 
% 
% % Set a weight for the manipulated variable rate
% mpcobj.Weights.ManipulatedVariablesRate = 1;

% Specify a lower bound for the manipulated variable.
mpcobj.ManipulatedVariables.Min = 0;

% The default output disturbance model contains a discrete-time integrator with dimensionless unity gain.
% Its input, |w|, is white noise with zero mean and unit variance. 
% Get the default disturbance model.
plant_output_disturbance = getoutdist(mpcobj);

% Augment the plant model with the disturbance model.
A = blkdiag(plant_discrete.A, plant_output_disturbance.A);
Bu = [plant_discrete.B; 0];
Cm = [plant_discrete.C plant_output_disturbance.C];
D = plant_discrete.D;

% Get the plant prediction model.
plant_prediction = ss(A, Bu, Cm, D, par.sample_time);

% Get the default measurement noise model for the output |y|, 
% which is a white noise with zero mean and unit variance.
plant_measurement_noise = ss(1, Ts = par.sample_time);

% Estimate |B| and |D| matrices.
B_est = [[plant_discrete.B; 0] [0; 0; plant_output_disturbance.B] [0; 0; 0]];
D_est = [plant_discrete.D plant_output_disturbance.D plant_measurement_noise.D];

% Obtain the noise covariance matrices |Q|, |R|, and |N|.
Q = [1.417059559107240e-10,-1.157557114736846e-07,0;-1.157557114736846e-07,9.455766804339990e-05,0;0,0,0.040000000000000];%B_est * B_est';
% Q = eye(3) * 1e-4;
R = D_est * D_est';
N = B_est * D_est';

% Obtain the gains, |L| and |M|, from the Kalman filter design.
G = eye(3);
H = zeros(1, 3);
[~, L, ~, M] = kalman(ss(A, [Bu G], Cm, [D H], par.sample_time), Q, R, N);

% Create a custom state estimation.
setEstimator(mpcobj, "custom")

%% Check the controllability of the model.
% Compute eigenvalues for stability analysis.
lambda = eig(A);

% Compute the controllability matrix.
controllability_matrix = ctrb(A, Bu); 
% Check the rank of the controllability matrix.
controllability_rank = rank(controllability_matrix);   

%% Open and simulate the model.
model = 'MPC_PMM_ml_gaslift';
open_system(model)
sim(model)

%% Plot figures.
fig = figure;

subplot(1, 2, 1)
plot(logsout{18}.Values.Time, logsout{18}.Values.Data, '--g', 'LineWidth', 1) % Injected gas flow rate (cm3/s)
xlabel('Time (s)')
ylabel('Q_{g, inj} (cm^3/s) (MV)', 'Interpreter', 'tex')
legend('MPC', 'Location', 'northoutside')
legend('boxoff')
title('')

subplot(1, 2, 2)
plot(logsout{27}.Values.Time, logsout{27}.Values.Data, '-k', 'LineWidth', 1) % Reference for liquid outflow rate (g/s)
hold on
plot(logsout{17}.Values.Time, logsout{17}.Values.Data, '--g', 'LineWidth', 1) % Liquid outflow rate (g/s)
xlabel('Time (s)')
ylabel('w_{l, out} (g/s) (OV)', 'Interpreter', 'tex')
legend('MPC',  'Location', 'northoutside', 'Orientation','horizontal')
legend('boxoff')
title('')

% savefig(fig, 'fig6_MPC_PMM_ml_Ts0_1_ph10_ch2_Q_Ts0_2_gaslift.fig') 
% close(fig)



