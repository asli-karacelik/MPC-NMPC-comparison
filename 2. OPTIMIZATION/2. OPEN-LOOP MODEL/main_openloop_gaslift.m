clear all
clc
%% Run the simulink file to obtain gas lift system's open-loop response

%% Parameters
par = function_parameters_openloop_gaslift;

%% Run the simulink file
plant_model = 'openloop_gaslift';
open(plant_model)
sim(plant_model)

%% Plot figure
f = figure;
plot(logsout{17}.Values.Time, logsout{17}.Values.Data, '-.b', 'LineWidth', 1) % Injected gas flow rate (m3/s)
legend('Q_{g, inj} = 2108.813609 cm^3/s', 'Location', 'northoutside')
legend('boxoff')
xlabel('Time (s)')
ylabel('w_{l,out} (g/s) (OV)', 'Interpreter', 'tex')
title('')
savefig(f, 'fig1_openloop_Qginjmax2108.813609_Ts0_05_gaslift.fig') 
close(f)

