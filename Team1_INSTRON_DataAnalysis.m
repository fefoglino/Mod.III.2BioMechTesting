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
for k = 1:length(materials)
    for n = 1:length(teamnames)
        data.(materials{k}).L0(1,n)   = (rawdata.(materials{k})(3,4*n-1));
        data.(materials{k}).Width(1,n) = (rawdata.(materials{k})(4,4*n-1));
        data.(materials{k}).Thickness(1,n)  = (rawdata.(materials{k})(5,4*n-1));
        data.(materials{k}).Area(1,n)   = (rawdata.(materials{k})(6,4*n-1));
        
        data.(materials{k}).Pos{1,n}    = (rawdata.(materials{k})(9:end,4*n-2));
        data.(materials{k}).Force{1,n}  = (rawdata.(materials{k})(9:end,4*n-1));
    end
end

%% Linearize data, calculate material properties
% I made this separate bc i kept getting weird index errors with the skeleton code
% and it was easier to debug like this
% basically any calcuations/modifications made to the original data end up
% in the linearized struct
for k = 1:length(materials)
    for n = 1:length(teamnames)
        % load the data for a specific material and team
        Force = cell2mat(data.(materials{k}).Force{1,n});
        Pos  = cell2mat(data.(materials{k}).Pos{1,n});

        % zero the data
        Force = Force - Force(1);
        Pos = Pos - Pos(1);

        % get rid of trailing NaN values only if they are present
        idx = find(isnan(Force), 1, 'first');
        if idx>0
            Force = Force(1:idx-1);
            Pos = Pos(1:idx-1);
        end

        % get rid of non-linear data only if max value comes before last value
        idx2 = find(Force == max(Force), 1, 'first');
        if idx2>0
            Force = Force(1:idx2-1);
            Pos = Pos(1:idx2-1);
        end

        % save the final versions of the force and pos arrays
        linearized.(materials{k}).Force{1,n} = Force;
        linearized.(materials{k}).Pos{1,n} = Pos;

        % find stiffness
        p = polyfit(Pos, Force, 1);
        linearized.(materials{k}).Stiffness(1,n) = p(1);

        % find stress
        A = cell2mat(data.(materials{k}).Area(1,n)); % note A is in mm2, aka 10^-6 m2
        Stress = Force./A; % so stress will be in MPa (10^6 Pa)
        linearized.(materials{k}).Stress{1,n} = Stress; 

        % find strain
        L0 = cell2mat(data.(materials{k}).L0(1,n));
        Strain = Pos./L0;
        linearized.(materials{k}).Strain{1,n} = Strain;

        % find extensibility
        linearized.(materials{k}).Extensibility(1,n) = 100.*max(Strain); % *100 to get a percentage

        % ultimate tensile strength
        linearized.(materials{k}).UltTenStrength(1,n) = max(Stress);

        % Young's modulus
        p2 = polyfit(Strain, Stress, 1);
        linearized.(materials{k}).E(1,n) = p2(1);

    end
    [stdStiff, meanStiff] = std(linearized.(materials{k}).Stiffness(1:end));
    [stdExt, meanExt] = std(linearized.(materials{k}).Extensibility(1:end));
    [stdUTS, meanUTS] = std(linearized.(materials{k}).UltTenStrength(1:end));
    [stdE, meanE] = std(linearized.(materials{k}).E(1:end));
end

%% Display Data:

%% Position and Force graphs
figure(1);
hold on
for k = 1:length(materials)
    subplot(2,2,k);

    for n = 1:length(teamnames)
        Force = linearized.(materials{k}).Force{1,n};
        Pos  = linearized.(materials{k}).Pos{1,n};

        plot(Pos, Force);
        hold on
        title(strcat("Force vs. Position (", materials{k}, ")"), 'Interpreter', 'None');
        xlabel('Position (mm)');
        ylabel('Force (N)');
    end
    legend(teamnames, 'Location', 'Best');
end
hold off

%% Stress and Strain graphs
figure(2);
hold on
for k = 1:length(materials)
    subplot(2,2,k);

    for n = 1:length(teamnames)
        Stress = linearized.(materials{k}).Stress{1,n};
        Strain = linearized.(materials{k}).Strain{1,n};

        plot(Strain, Stress);
        hold on
        title(strcat("Stress vs. Strain (", materials{k},")"), 'Interpreter', 'None');
        xlabel("Strain");
        ylabel("Stress (MPa)");
    end
    legend(teamnames, 'Location', 'Best');
end

%% Statistical analysis

figure(3);
hold on
for k = 1:length(materials)
    subplot(2,2,k);
    labels = {"Young's Modulus", "Stiffness", "Extensibility", "Ultimate Tensile Strength"};
    boxMatrix = cat(1, linearized.(materials{k}).E(1:end), linearized.(materials{k}).Stiffness(1:end), linearized.(materials{k}).Extensibility(1:end), linearized.(materials{k}).UltTenStrength(1:end));
    boxchart(boxMatrix);
end

%% Useful Commands:
    % nanmean(cell2mat(data.Skin_Raw.Height)) --> Calculates the mean when
    % there are NaN entries in a matrix

    % gradient() --> computes the derivative of a dataset