%% Run and plot scenarios
function runPlot_SA_scenarios(model,data,args)
% runPlot_SA_scenarios  Simulate and plot a one-at-a-time parameter sweep.
%   runPlot_SA_scenarios(model,data,args) builds a scenario that sweeps the
%   parameter args.ParaVar over args.ParaVarValues (group 111111111), simulates
%   the model, and plots the CSF protein-percent response for each swept value
%   (overlaying matching data).
%
%   Inputs:
%     model - SimBiology model.
%     data  - dataset table for data overlays.
%     args  - struct of settings (ParaVar, ParaVarValues, figure options).

%% User input for the runPlot
stopTime = hours(days(320)); % simulation duration
alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations

%% set up the scenarios
group=[111111111];
% color_map = ["blue","red","yellow","purpledark","green","cyan","reddark","teal","brown","pink"];
color_map = ["purple1","purple2","purple3","purple4","black"];
p = sbioselect(model,'Name','group','type','parameter');
if isempty(p); p = addparameter(model,'group',0,'Units','dimensionless'); end

% get tasks with the same order as group
[Lia,Locb] = ismember(args.taskTable.Group,group);
tasks = args.taskTable(Lia,:);
[~,p] = sort(Locb(Lia));
tasks = tasks(p,:);

% Display Name from the task spreadsheet
displayName = tasks.DisplayName;

% Add group similar to gQSPsim task setups
sObj = SimBiology.Scenarios;
add(sObj,'elementwise','group',tasks.Group);

% Add Variants as elementwise
varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
for i=varCol
    add(sObj,'elementwise','variants',tasks.(tasks.Properties.VariableNames{i})(:));
end

% Add doses as elementwise
doseCol = find(contains(lower(tasks.Properties.VariableNames), lower('Dose')));
for i=doseCol
    add(sObj,'elementwise','dose',tasks.(tasks.Properties.VariableNames{i})(:));
end

% Add cartesian variability
if args.variability
    sObj = add_eta_BOAT_v0(model,sObj,args);
end

% Add SA parameter
add(sObj,'Cartesian',args.ParaVar,args.ParaVarValues);


% Build and examine the scenario table
scenariosTbl = generate(sObj);
sObj.RandomSeed = sum(double('jafarnem'));

%% Get list of states to log for plotting
% For scenarios where sim time efficiency matters, we can just log the
% needed states and reduce memory usage and potentially increase speed

cs = getconfigset(model);
states = cs.RuntimeOptions.StatesToLog;
for i = 1:length(states)
    observableNames{i} = [states(i).PartiallyQualifiedName];
end

%% Run the model using the simFunc
tic % start the timer

simFunc = createSimFunction(model,sObj,observableNames,[],[],'UseParallel',args.UseParallel, 'AutoAccelerate',1);

simDataArray = simFunc(sObj, stopTime, [], 0:1:stopTime);

elapsedTime = toc;  % stop the timer and calculate elapsed time
fprintf('The simulation time is %.2f seconds.\n', elapsedTime);

%% Add any post simulation observables that need the timeseries like min/max
tic

% simDataArray = addobservable(simDataArray, 'Cmax_TDB_Cen_ugmL', 'max(TDB_Cen_ugmL)');

elapsedTime = toc;  % stop the timer and calculate elapsed time
fprintf('The observable calculation time is %.2f seconds.\n', elapsedTime);

%% Plot PK
% create plots 


f = figure;

for i=1:length(group)
    
    for ii=1:length(args.ParaVarValues)
         
        f.Position = args.figSize(1,:);
        colororder("gem12");
        colorOrig = get(gca,'colororder');
        color = colorOrig;
    % Total Oligo Plasma
    sp = subplot(1,1,1); hold on
    sp.Box="on";

    % plot the simulation
    speciesModel = 'V_csf.Prot_perc';
    simDataTemp = simDataArray(ii);
    simDataTemp = selectbyname(simDataTemp,speciesModel);
    ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_map(ii)), DisplayName=args.ParaVarValues(ii));
    
    dataTemp = getData(data, Group=250813004, Analyte='Protein', Matrix={'CSF'});
    if ~isempty(dataTemp)
        % plotData: gets the filtered time vs conc data with option for
        % different types of plots
        % TASK: make a plot XY version of the plotData function
        % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
        ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker(color_map(3)), LineStyle="none");
        
    end
    % make the plot look nice
    ax.Title.String = displayName{i};
    ax.FontSize = args.fontSize;
    % y-Axis
    ax.YLim = [0, 110];
    yticks(ax,[0 25 50 75 100]);
    % ax.YScale = 'log';
    ax.YLabel.String = "Protein %";
    % x-Axis
    xTimeAxis(ax,XTimeUnit="day")
    xticks(ax,[0:1440:7200]);
    xticklabels(ax,{'0','60','120','180','240','300'});
    % yAxisSecondUnit(ax,unitLeft='nM',unitRight='ug/mL',MW=150);

    hold on;

    % f = figure(2); 
    % f.Position = args.figSize(1,:);
    % colororder("gem12");
    % colorOrig = get(gca,'colororder');
    % color = colorOrig;
    % 
    % % Total Ab Plasma
    % sp = subplot(1,1,1); hold on
    % sp.Box="on";
    % 
    % % plot the simulation
    % speciesModel = 'body.totalAntibody_plasma_nM';
    % simDataTemp = simDataArray(ii);
    % simDataTemp = selectbyname(simDataTemp,speciesModel);
    % ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_map(ii)), DisplayName=args.ParaVarValues(ii));
    % 
    % dataTemp = getData(data, Group=250014807, Analyte='PK', Matrix={'Plasma'});
    % if ~isempty(dataTemp)
    %     % plotData: gets the filtered time vs conc data with option for
    %     % different types of plots
    %     % TASK: make a plot XY version of the plotData function
    %     % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
    %     ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('blue'), LineStyle="none");
    % end
    % 
    % % make the plot look nice
    % ax.Title.String = displayName{i};
    % ax.FontSize = args.fontSize;
    % % y-Axis
    % ax.YLim = [1e-2, 5e3];
    % ax.YScale = 'log';
    % ax.YLabel.String = "Total Ab (nM)";
    % % x-Axis
    % xTimeAxis(ax,XTimeUnit="day")
    % % yAxisSecondUnit(ax,unitLeft='nM',unitRight='ug/mL',MW=150);
    % 
    % hold on;
    % 
    % f = figure(3); 
    % f.Position = args.figSize(1,:);
    % colororder("gem12");
    % colorOrig = get(gca,'colororder');
    % color = colorOrig;
    % 
    % % Total Oligo Liver
    % sp = subplot(1,1,1); hold on
    % sp.Box="on";
    % 
    % % plot the simulation
    % speciesModel = 'body.totalASO_liver_tot_nM';
    % simDataTemp = simDataArray(ii);
    % simDataTemp = selectbyname(simDataTemp,speciesModel);
    % ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_map(ii)), DisplayName=args.ParaVarValues(ii));
    % 
    % dataTemp = getData(data, Group=250014807, Analyte='OligoPK', Matrix={'Liver'});
    % if ~isempty(dataTemp)
    %     % plotData: gets the filtered time vs conc data with option for
    %     % different types of plots
    %     % TASK: make a plot XY version of the plotData function
    %     % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
    %     ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('red'), LineStyle="none");
    % end
    % % make the plot look nice
    % ax.Title.String = displayName{i};
    % ax.FontSize = args.fontSize;
    % % y-Axis
    % ax.YLim = [1e2, 1e4];
    % % ax.YScale = 'log';
    % ax.YLabel.String = "Total Oligo (nM)";
    % % x-Axis
    % xTimeAxis(ax,XTimeUnit="day")
    % % yAxisSecondUnit(ax,unitLeft='nM',unitRight='ug/mL',MW=150);
    % 
    % hold on;

    % f = figure(4); 
    % f.Position = args.figSize(1,:);
    % colororder("gem12");
    % colorOrig = get(gca,'colororder');
    % color = colorOrig;
    % 
    % % Free Oligo Brain
    % sp = subplot(1,1,1); hold on
    % sp.Box="on";
    % 
    % % plot the simulation
    % speciesModel = 'body.freeASO_brain_nM';
    % simDataTemp = simDataArray(ii);
    % simDataTemp = selectbyname(simDataTemp,speciesModel);
    % ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_map(ii)), DisplayName=args.ParaVarValues(ii));
    % 
    % % make the plot look nice
    % ax.Title.String = displayName{i};
    % ax.FontSize = args.fontSize;
    % % y-Axis
    % ax.YLim = [1e-2, 1e2];
    % ax.YScale = 'log';
    % ax.YLabel.String = "Free Oligo (nM)";
    % % x-Axis
    % xTimeAxis(ax,XTimeUnit="day")
    % % yAxisSecondUnit(ax,unitLeft='nM',unitRight='ug/mL',MW=150);
    % 
    % hold on;
    end
end


%% saving the whole figure here
% for some builds it makes sense to save individual axis separately and
% also build on them. exportgraphics can take ax as input and the naming
% can be adjustes to include specific treatment, dose, matrix or analyte
% infomration of interest
% if want to add date: [char(datetime("today"),'yyyyMMdd_'),dbstack().name,'_']
% funcName = char(dbstack().name);
% append(dbstack().name,'_test',args.figFormat)
