clear all
clc
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

%% Get parameters
par = function_parameters_MPC_MD_gaslift;
% % Store parameters in an Excel file.
% T = {'Run', 'Sample Time (s)', 'Reservoir Pressure Change [g/s])', 'Comment'};
% writecell(T, 'parameters_MPC_MD_gaslift.xlsx', 'Sheet', 1, 'Range', 'A1:D1');
% T = table({'1'}, par.sample_time, par.wloutref);
% writetable(T, 'parameters_MPC_MD_gaslift.xlsx', 'Sheet', 1, 'WriteVariableNames', 0, 'Range', 'A2:C2');
% winopen('parameters_MPC_MD_gaslift.xlsx')

%% Linearization
% The internal plant model is obtained by linearizing a nonlinear plant model,
% |model_MD_gaslift|, in Simulink at a nonzero steady-state operating point.

% Open and simulate the nonlinear model. 
plant_model = 'model_MD_gaslift';
sim(plant_model)

% Generate an operating point specification for the model initial conditions.
operating_point_specification = operspec(plant_model);
% Define the measured disturbance as a known variable.
operating_point_specification.Inputs(2).Known = 1;
% Set an initial condition for the measured disturbance.
operating_point_specification.Inputs(2).u = par.pres; 

% Obtain the operating point from the model's initial conditions.
[operating_point, operating_report] = findop(plant_model, operating_point_specification);

% Specify nominal state, output, and input values from the operating point.
x0 = [operating_report.States(1).x; operating_report.States(2).x];
y0 = operating_report.Outputs.y;
u0 = [operating_report.Inputs(1).u; par.pres];
          
% Linearize the plant at the operating point.
plant = linearize(plant_model, operating_point ); 
% Introduce the manipulated variable and the measured disturbance into the linearized model.
plant = setmpcsignals(plant, MV = 1, MD = 2);

%% Configure MPC 
% mpcobj = mpc(plant, sample time, prediction horizon steps, control horizon steps, weights, manipulated variables, output variables, disturbance variables)
% mpcobj = mpc(plant, ts, P, M, W, MV, OV, DV)
mpcobj = mpc(plant, par.sample_time, par.prediction_horizon, par.control_horizon);

% Set the nominal values in the controller.
mpcobj.Model.Nominal = struct('X', x0, 'U', u0, 'Y', y0);

% Change the plant model from continuous to discrete time.
plant_discrete = c2d(plant, par.sample_time);

% % Set minimum and maximum manipulated variable rates (injected gas flow rate [m3/s]).
% mpcobj.ManipulatedVariables = struct('RateMin', -100, 'RateMax', 100); 
% % Set a weight for the manipulated variable rate
% mpcobj.Weights.ManipulatedVariablesRate = 1;

% Specify a lower bound for the manipulated variable.
mpcobj.ManipulatedVariables(1).Min = 0;

% % Scaling
% mpcobj.DisturbanceVariables.ScaleFactor = 1e7;

% The default output disturbance model contains a discrete-time integrator with dimensionless unity gain.
% Its input, |w|, is white noise with zero mean and unit variance. 
% Get the default disturbance model.
plant_output_disturbance = getoutdist(mpcobj);

% Augment the plant model with the disturbance model.
A = blkdiag(plant_discrete.A, plant_output_disturbance.A);
Bu = [plant_discrete.B; 0 0];
Cm = [plant_discrete.C, plant_output_disturbance.C];
D = plant_discrete.D;

% Get the plant prediction model.
plant_prediction = ss(A, Bu, Cm, D, par.sample_time);

% Get the default measurement noise model for the output |y|, 
% which is a white noise with zero mean and unit variance.
plant_measurement_noise = ss(1, Ts = par.sample_time);

% Estimate |B| and |D| matrices.
B_est = [[plant_discrete.B; 0 0] [0; 0; plant_output_disturbance.B] [0; 0; 0]];
D_est = [plant_discrete.D, plant_output_disturbance.D, plant_measurement_noise.D];

% Obtain the noise covariance matrices |Q|, |R|, and |N|.
Q = B_est * B_est';
R = D_est * D_est';
N = B_est * D_est';

% Obtain the gains, |L| and |M|, from the Kalman filter design.
G = eye(3);
H = zeros(1, 3);
% [kalmf,L,P,Mx,Z,My] = kalman(___) 
% kalmf — Steady-state filter model
% P(continuous), Z(discrete) — Steady-state error covariances
[~, L, ~, M] = kalman(ss(A, [Bu G], Cm, [D H], par.sample_time), Q, R, N); % G and H is for disturbance

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
model = 'MPC_MD_gaslift';
open_system(model)
sim(model)

%% Plot figures.
fig = figure;
subplot(1, 2, 1)
yyaxis right
plot(logsout{18}.Values.Time, logsout{18}.Values.Data, '--g', 'LineWidth', 1) % Injected gas flow rate (cm3/s)
xlabel('Time (s)')
ylabel('Q_{g, inj} (cm^3/s) (MV)', 'Interpreter', 'tex')
% legend('MPC', 'Location', 'northoutside')
% legend('boxoff')
title('')

subplot(1, 2, 2)
plot(logsout{27}.Values.Time, logsout{27}.Values.Data, '-k', 'LineWidth', 1) % Reference for liquid outflow rate (g/s)
hold on
plot(logsout{17}.Values.Time, logsout{17}.Values.Data, '--g', 'LineWidth', 1) % Liquid outflow rate (g/s)
xlabel('Time (s)')
ylabel('w_{l, out} (g/s) (OV)', 'Interpreter', 'tex')
legend('Reference', 'MPC', 'Location', 'northoutside', 'Orientation','horizontal')
legend('boxoff')
title('')

% savefig(fig, 'fig6_MPC_MD_0_05pres_Ts0_1_ph10_ch2_gaslift.fig') 
% close(fig)
