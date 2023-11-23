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
end

%% Display Data:

% this makes the graph titles/labels look better
matLabels = {'Silicone, Thin', 'Silicone, Thick', ...
             'Skin, Raw', 'Skin, Treated'};

%% Position and Force graphs
figure(1);
hold on
for k = 1:length(materials)
    subplot(2,2,k);

    for n = 1:length(teamnames)
        Force = linearized.(materials{k}).Force{1,n};
        Pos  = linearized.(materials{k}).Pos{1,n};

        plot(Pos, Force, 'LineWidth', 1);
        hold on
        title(strcat("Force vs. Position (", matLabels{k}, ")"), 'Interpreter', 'None');
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

        plot(Strain, Stress, 'LineWidth', 1);
        hold on
        title(strcat("Stress vs. Strain (", matLabels{k},")"), 'Interpreter', 'None');
        xlabel("Strain");
        ylabel("Stress (MPa)");
    end
    legend(teamnames, 'Location', 'Best');
end

%% Statistical analysis

%% find mean, std for properties and init measurements
for k = 1:length(materials)

    % material properties
    [stdStiff, meanStiff] = std(linearized.(materials{k}).Stiffness(1:end));
    [stdExt, meanExt] = std(linearized.(materials{k}).Extensibility(1:end));
    [stdUTS, meanUTS] = std(linearized.(materials{k}).UltTenStrength(1:end));
    [stdE, meanE] = std(linearized.(materials{k}).E(1:end));
    linearized.(materials{k}).YM_St_Ex_UTS(1:4) = [meanE meanStiff meanExt meanUTS];
    linearized.(materials{k}).STD_YM_St_Ex_UTS(1:4) = [stdE stdStiff stdExt stdUTS];

    % thickness, width, initial length
    [stdThick, meanThick] = std(cell2mat(data.(materials{k}).Thickness(1:end)));
    [stdWidth, meanWidth] = std(cell2mat(data.(materials{k}).Width(1:end)));
    [stdL0, meanL0] = std(cell2mat(data.(materials{k}).L0(1:end)));
    linearized.(materials{k}).ThWdL0(1:3) = [meanThick meanWidth meanL0];
    linearized.(materials{k}).STD_ThWdL0(1:3) = [stdThick stdWidth stdL0];

end

%% ANOVA

% empty matrices to compile data in correct format for ANOVA1()
EMat = [];
StiffMat = [];
ExtMat = [];
UTSMat = [];

% put data into matrices in correct format
for k = 1:length(materials)
    EMat = [EMat, linearized.(materials{k}).E(1:end)'];
    StiffMat = [StiffMat, linearized.(materials{k}).Stiffness(1:end)'];
    ExtMat = [ExtMat, linearized.(materials{k}).Extensibility(1:end)'];
    UTSMat = [UTSMat, linearized.(materials{k}).UltTenStrength(1:end)'];
end

% perform ANOVA
[pE, tblE, statE] = anova1(EMat, materials, 'off');
[pStiff, tblStiff, statStiff] = anova1(StiffMat, materials, 'off');
[pExt, tblExt, statExt] = anova1(ExtMat, materials, 'off');
[pUTS, tblUTS, statUTS] = anova1(UTSMat, materials, 'off');

% multiple comparison tests
[compE, mE, hE, gE] = multcompare(statE, 'Display', 'off');
[compStiff, mStiff, hStiff, gStiff] = multcompare(statStiff, 'Display', 'off');
[compExt, mExt, hExt, gExt] = multcompare(statExt, 'Display', 'off');
[compUTS, mUTS, hUTS, gUTS] = multcompare(statUTS, 'Display', 'off');

% comp = [group1, group2, CI1, dif in mean, CI2, p]

compAll = {compE, compStiff, compExt, compUTS};

%% Material properties box plots

labels = {"Young's Modulus", "Stiffness", "Extensibility", "Ultimate Tensile Strength (MPa)"};
for i = 1:length(labels)
    figure(i+3);
    hold on
    barMatrix = [];
    errorUpper = [];

    for k = 1:length(materials)
        % make a 1x4 matrix with the means of the values for each material
        % in order
        oldMatrix = barMatrix;
        newData = linearized.(materials{k}).YM_St_Ex_UTS(i);
        barMatrix = [oldMatrix newData];
        % 1x4 matrix of the STdev of the values for each material
        oldErr = errorUpper;
        newErr = linearized.(materials{k}).STD_YM_St_Ex_UTS(i);
        errorUpper = [oldErr newErr];
    end

    % set the minimum error bar value to 0
    errorLower = errorUpper;
    for n = 1:length(errorUpper)
        if barMatrix(n)-errorLower(n)<0
            errorLower(n) = barMatrix(n);
        end
    end

    barGraph = bar(matLabels, barMatrix);
    title(labels{i});

    % add mean and std labels
    for m = 1:length(labels)
        toDisplay = {sprintf("%.4g", barMatrix(m)), " \pm " + sprintf("%.4g", errorUpper(m))};
        text(barGraph.XEndPoints(m), barGraph.YEndPoints(m), toDisplay, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
    end
    
    % add color-coded brackets for p vals btwn each group
    % this part took me longer than the rest combined
    bct = 0.9;
    colors = {"#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE"};
    for q = 1:length(compAll{i})
        if compAll{i}(q, 6)<0.05 % add bracket if p<0.05 for a particular comparison
            % graphDim = [0.13 0.11 0.775 0.815]
            g1 = 0.2.*compAll{i}(q,1);
            g2 = 0.2.*compAll{i}(q,2);
            xl1 = 0.775.*g1 + 0.13;
            xl2 = 0.775.*g2 + 0.13;
            % three lines = one bracket
            annotation('line', [xl1 xl2], [bct bct], 'LineWidth', 2, 'Color', colors{q});
            annotation('line', [xl1 xl1], [bct 0.8], 'LineWidth', 2, 'Color', colors{q});
            annotation('line', [xl2 xl2], [bct 0.8], 'LineWidth', 2, 'Color', colors{q});
            pVal = sprintf('p = %.4f', compAll{i}(q,6));
            annotation('textbox', [0.1 bct 0.1 0.1], 'String', pVal, 'Color', colors{q}, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 16);
            bct = bct - 0.02;
        end
    end

    % add error bars
    hold on
    stdMarkers = errorbar(1:length(materials), barMatrix, errorLower, errorUpper);
    stdMarkers.Color = [0 0 0];
    stdMarkers.LineStyle = 'none';
    hold off
end

