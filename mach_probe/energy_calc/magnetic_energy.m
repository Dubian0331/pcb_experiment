%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters and Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
u0 = 4 * pi * 10^(-7); % Permeability of free space

% File paths
% caseI
% date = 240610;
% startFile = 4289;
% endFile = 4319;

% caseO
date = 240611;
startFile = 4363;
endFile = 4389;

dataFolder = 'G:\My Drive\lab\lab_data\pre_processed\';
saveFolder = 'G:\My Drive\lab\lab_data\mach_probe\energy_calc\';
saveFolderPath = strcat(saveFolder, '\', num2str(date));



% User-defined time range (supports scalar list or range)
trange_user = 470:520; % Specify the time points or range you want to analyze

% Variables for storing results
total_B_all = []; % To store total B values for all files
trange_actual = []; % To store the actual trange values corresponding to user input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for fileNum = startFile:endFile
    % Construct file name and check existence
    fileName = sprintf('a039_%d.mat', fileNum);
    fullPath = fullfile(dataFolder, fileName);
    if exist(fullPath, 'file') == 2
        % Load data
        load(fullPath, 'grid2D', 'data2D');
        
        % Retrieve trange from the data and reshape to a vector
        trange_data = data2D.trange(:); % Ensure it's a column vector
        
        % Find indices corresponding to trange_user
        idx_user = find(ismember(trange_data, trange_user)); % Match user-defined range
        trange_actual = trange_data(idx_user); % Store actual matching values
        
        if isempty(idx_user)
            warning('No valid time points found in trange for file: %s', fileName);
            continue;
        end
        
        % Initialize total_B for current file
        total_B = zeros(1, length(idx_user));
        
        % Loop through the matched time points
        for i = 1:length(idx_user)
            t_idx = idx_user(i);
            
            % Extract magnetic field data
            Bz = data2D.Bz(:, :, t_idx); % Magnetic field in z direction
            Br = data2D.Br(:, :, t_idx); % Magnetic field in r direction
            Bt = data2D.Bt(:, :, t_idx); % Magnetic field in t direction
            
            % Calculate the total magnetic field energy density
            B_total = sum(Bz.^2 + Br.^2 + Bt.^2); % Sum over spatial range
            total_B(i) =  B_total / (2 * u0); % Magnetic energy density
        end
        
        % Append results for plotting
        total_B_all = [total_B_all; total_B]; % Add total_B for this file
    else
        warning('File not found: %s', fullPath);
    end
end

% Calculate mean and standard deviation over all files
mean_total_B = mean(total_B_all, 1); % Average across files
std_total_B = std(total_B_all, 0, 1); % Standard deviation across files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Plot mean with error bars
e = errorbar(trange_actual, mean_total_B, std_total_B, '-o');
xlabel('Time [\mus]', 'FontSize', 20);
ylabel('Total Magnetic Energy Density [J]', 'FontSize', 20);
title('Time Evolution of Total Magnetic Energy Density', 'FontSize', 24);

ax = gca;
ax.FontSize = 16;
e.LineWidth = 2;
e.MarkerSize = 8;

% Save plot
if ~exist(saveFolderPath, 'dir')
    mkdir(saveFolderPath);
end
saveas(gcf, fullfile(saveFolderPath, sprintf('%d_%d_%s_time_magnetic_energy.png', startFile, endFile, num2str(date))));

disp('Time evolution of total magnetic energy density plot completed.');


%%
total_B_energy = sum(total_B_all(:))



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot at the same time. caseI and caseO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
u0 = 4 * pi * 10^(-7); % Permeability of free space

% Case configurations
cases = {struct('name', 'CaseI', 'shots', 4289:4319, 'color', 'r'), ...
         struct('name', 'CaseO', 'shots', 4363:4389, 'color', 'b')};

% File paths
dataFolder = 'G:\My Drive\lab\lab_data\pre_processed\';
saveFolder = 'G:\My Drive\lab\lab_data\mach_probe\energy_calc\';
saveFolderPath = strcat(saveFolder, '\', 'Comparison_Case-I_Case-O');

% Time range
trange_user = 470:520; % Specify the time points or range you want to analyze

% Initialize variables for results
results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Processing for Each Case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for c = 1:length(cases)
    caseData = cases{c};
    shots = caseData.shots;
    total_B_all = [];
    trange_actual = [];
    
    % Process each file in the case
    for shot = shots
        % Construct file name and check existence
        fileName = sprintf('a039_%d.mat', shot);
        fullPath = fullfile(dataFolder, fileName);
        if exist(fullPath, 'file') == 2
            % Load data
            load(fullPath, 'grid2D', 'data2D');
            
            % Retrieve trange from the data and reshape to a vector
            trange_data = data2D.trange(:); % Ensure it's a column vector
            
            % Find indices corresponding to trange_user
            idx_user = find(ismember(trange_data, trange_user)); % Match user-defined range
            trange_actual = trange_data(idx_user); % Store actual matching values
            
            if isempty(idx_user)
                warning('No valid time points found in trange for file: %s', fileName);
                continue;
            end
            
            % Initialize total_B for current file
            total_B = zeros(1, length(idx_user));
            
            % Loop through the matched time points
            for i = 1:length(idx_user)
                t_idx = idx_user(i);
                
                % Extract magnetic field data
                Bz = data2D.Bz(:, :, t_idx); % Magnetic field in z direction
                Br = data2D.Br(:, :, t_idx); % Magnetic field in r direction
                Bt = data2D.Bt(:, :, t_idx); % Magnetic field in t direction
                
                % Calculate the total magnetic field energy density
                B_total_squared = sum(sum(Bz.^2 + Br.^2 + Bt.^2)); % Sum over spatial range
                total_B(i) = B_total_squared / (2 * u0); % Magnetic energy density
            end
            
            % Append results for plotting
            total_B_all = [total_B_all; total_B]; % Add total_B for this file
        else
            warning('File not found: %s', fullPath);
        end
    end
    
    % Calculate mean and standard deviation over all files
    mean_total_B = mean(total_B_all, 1); % Average across files
    std_total_B = std(total_B_all, 0, 1); % Standard deviation across files
    
    % Store results
    results = [results, struct('name', caseData.name, 'color', caseData.color, ...
                               'trange', trange_actual, 'mean', mean_total_B, ...
                               'std', std_total_B)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on;

% Plot data for each case
for c = 1:length(results)
    result = results(c);
    errorbar(result.trange, result.mean, result.std, '-o', ...
        'Color', result.color, 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', result.name);
end

xlabel('Time [\mus]', 'FontSize', 20);
ylabel('Total Magnetic Energy Density [J]', 'FontSize', 20);
ylim([0 9e5]);
title('Comparison of Case-I and Case-O: Time Evolution of Magnetic Energy', 'FontSize', 24);
legend('Location', 'best');
grid on;

ax = gca;
ax.FontSize = 16;

% Save plot
if ~exist(saveFolderPath, 'dir')
    mkdir(saveFolderPath);
end
saveas(gcf, fullfile(saveFolderPath, 'Comparison_Case-I_Case-O_Magnetic_Energy.png'));

disp('Time evolution of total magnetic energy density comparison completed.');
