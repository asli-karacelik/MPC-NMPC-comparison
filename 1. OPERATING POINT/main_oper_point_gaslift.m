clear all
clc
%% Parameters
par = function_parameter_oper_point_gaslift;

%% Operating point specification for the current model initial condition
plant_model = 'model_oper_point_gaslift';
open_system(plant_model)
sim(plant_model)
operating_point_specification = operspec(plant_model);

%% The operating point for this initial condition
[operating_point, operating_report] = findop(plant_model, operating_point_specification);

%% Nominal state, output, and input values from the computed operating point
x0 = [operating_report.States(1).x; operating_report.States(2).x];
y0 = operating_report.Outputs.y;
u0 = operating_report.Inputs.u;