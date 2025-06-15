clear all
clc
tic
%% Description: 
    % A vertical pipe is submerged in a water-filled reservoir. When the air is 
    % injected into the bottom of the pipe, water starts to flow out. However,  
    % without air injection, there is no water outflow from the pipe. 

% Input
    % Qginj: injected gas flow rate (cm3/s) [manipulated variable]
    
% States 
    % mg: gass mass in the pipe (g)
    % ml: liquid mass in the pipe (g)

% Output 
    % wlout: liquid outflow rate (g/s)

%% Folder and file settings
folder = '.'; % Current directory
file_pattern = '*.csv'; % Match the file extension
file_list = dir(fullfile(folder, file_pattern)); % Get list of files

% Excel file loop
for file_idx = 1:length(file_list)
    % Get the filename
    file_name = file_list(file_idx).name;
    
    % Extract submergence ratio from the filename
    S_str = regexp(file_name, '(\d+\_\d+)', 'tokens'); 

    % Extract the first matched token
    S_str = S_str{1}{1};  % 'tokens' returns a cell of cells

    % Convert the string with underscore to a numeric value
    S_str = strrep(S_str, '_', '.');  % Replace underscores with a dot
    S = str2double(S_str);  % Convert the string to a numeric value

    % Load parameters with dynamic submergence ratio (S)
    par = function_parameter_Stenning(S);
    
    % Read data from the current CSV file
    file_path = fullfile(folder, file_name);
    data = readtable(file_path);
    air = data{:, 1} * 28316.85; % Convert air flow rates (ft3/s to cm3/s)
    water = data{:, 2} * (28316.85 * par.rhol); % Convert water flow rates (ft3/s to g/s)

    % Set bounds for air flow rates
    air_min = 0;
    air_max = 110;%max(air);
    increment = (air_max - air_min) / 100000;
    
    %% Fsolve computations
    results_wlout = [];
    results_exitflag = [];
    x1 = par.mg0; % Initial gas mass (g)
    x2 = par.ml0; % Initial liquid mass (g)
    x0 = [x1 x2];
    
    for Qginj = air_min:increment:air_max % Volumetric air inflow rate (cm3/s)
        opt = optimoptions('fsolve', 'MaxIter', 1e3, 'MaxFunEvals', 1e3, 'Display','off');
        [x, fval, exitflag, output] = fsolve(@(x)function_model_Stenning(x, Qginj, par), x0, opt);
        [F, wlout] = function_model_Stenning(x, Qginj, par);
        results_wlout = [results_wlout wlout];
          
    end
    %% Plot the results
    Qginj = air_min:increment:air_max;
    figure;
    plot(Qginj, results_wlout, '-.g', 'LineWidth', 1);
    hold on;
    plot(air, water, '-k');
    xlabel('Q_{g,inj} (cm^3/s)');
    ylabel('w_{l,out} (g/s)');
    legend('Theory', 'Experiment');
    title(['Results for Submergence Ratio = ', num2str(S)]);
    hold off;
end

%% Calculate optimal point
[max_wlout, idx] = max(results_wlout);
Qginj_max = Qginj(idx);
fprintf('Maximum w_{l,out} = %.6f g/s at Q_{g,inj} = %.6f cm^3/s\n', max_wlout, Qginj_max);
toc