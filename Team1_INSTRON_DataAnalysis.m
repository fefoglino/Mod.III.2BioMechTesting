%% BMEN E3810: Biomedical Engineering Laboratory I (INSTRON Data Analysis)
%  Written by: Dr. Lauren N. Heckelman
%  Fall 2023

%% Initialize the Workspace:
clear; clc; close all;

%% Define Data Parameters:
materials = {'Silicone_Thin', 'Silicone_Thick', ...
             'Skin_Raw', 'Skin_Treated'};
teamnames = {'Wed01',  'Wed02',  'Wed03',  'Wed04', ...
             'Wed05',  'Wed06', 'Wed07', ...
             'Thurs01','Thurs02','Thurs03', 'Thurs04', ...
             'Thurs05','Thurs06','Thurs07'};

filename = 'III.2 - Biomechanical Testing - Class Data.xlsx';

%% Import Raw Data:
for k = 1:length(materials)
    rawdata.(materials{k}) = readcell(filename, 'Sheet', materials{k});
    rawdata.(materials{k})(cellfun(@(x) any(ismissing(x)), rawdata.(materials{k}))) = {NaN};
end

%% Organize Data:
for k = 1:length(materials)-2
    for n = 1:length(teamnames)
        data.(materials{k}).L0(1,n)   = (rawdata.(materials{k})(3,4*n-1));
        data.(materials{k}).Width(1,n) = (rawdata.(materials{k})(4,4*n-1));
        data.(materials{k}).Thickness(1,n)  = (rawdata.(materials{k})(5,4*n-1));
        data.(materials{k}).Area(1,n)   = (rawdata.(materials{k})(6,4*n-1));
        
        data.(materials{k}).Pos{1,n}    = (rawdata.(materials{k})(9:end,4*n-2));
        data.(materials{k}).Force{1,n}  = (rawdata.(materials{k})(9:end,4*n-1));
    end
end

%% Display Data:
for k = 1:length(materials)
    figure(k);

    for n = 1:length(teamnames)
        Force = cell2mat(data.(materials{k}).Force{1,n});
        Pos  = cell2mat(data.(materials{k}).Pos{1,n});
        plot(Pos, Force);
        hold on
        title(strcat("Force_vs_Position_(", materials{k}, ")"), 'Interpreter', 'None');
        xlabel('Position (mm)');
        ylabel('Force (N)');
    end
    legend(teamnames, 'Location', 'Best');
end

%% Useful Commands:
    % nanmean(cell2mat(data.Skin_Raw.Height)) --> Calculates the mean when
    % there are NaN entries in a matrix

    % gradient() --> computes the derivative of a dataset

%% EXAMPLE OF SKIN_RAW PARTIAL ANALYSIS:
    figure(7); clf;

    for n = 1:length(teamnames)
        Force = cell2mat(data.(materials{3}).Force{1,n});
        Pos  = cell2mat(data.(materials{3}).Pos{1,n});
        
        % Invert force/position data 
        Force = -1.*Force;
        Pos = -1.*Pos;

        % Zero force/position data
        Force = Force - Force(1);
        Pos = Pos - Pos(1);
        
        % Determine if and where NaN entries exist in force/position data
        idx = find(isnan(Force), 1, 'first');
        Force = Force(1:idx-1);
        Pos = Pos(1:idx-1);
    
        % Crop off all values outside of linear range (You need to figure
        % this out! :))
        idx2 = find(Force == max(Force), 1, 'first');
        Force_linear = Force(:);
        Pos_linear = Pos(:);

        % Calculate linear regression of force/position data to determine
        % stiffness
        p = polyfit(Pos_linear, Force_linear, 1);
        Stiffness(1,n) = p(1);

        % Plot force/position data
        plot(Pos_linear, Force_linear);
        hold on
        title(strcat("Force_vs_Position_(", materials{3}, ")"), 'Interpreter', 'None');
        xlabel('Position (mm)');
        ylabel('Force (N)');
    end
    legend(teamnames, 'Location', 'Best');