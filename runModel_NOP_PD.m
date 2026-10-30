function runModel_NOP_PD(model,data,args)
% runModel_NOP_PD  PD scenario: MAPT, APP siRNA, HTT, and Jacifusen protein knockdown.
%   Builds and applies the scenario variants and doses, simulates the model
%   for each PD program (MAPT, APP siRNA, HTT groups, Jacifusen), and plots
%   the CSF/brain protein percent-baseline panels.
%
%   Inputs:
%     model - SimBiology model.
%     data  - dataset table for data overlays.
%     args  - struct of settings.
%% Plot PK
nCol = 4; % 4;
nRow = 2; % 3
nFig = 9; % ceil(length(group)/nCol);



% create plots 
f = figure(1); 
f.Position = args.figSize(1,:);
colororder("gem12");
colorOrig = get(gca,'colororder');
color = colorOrig;
% TASK: add color name to the task table and create a color selector that picks the rgb values based on the name
% modify the colororder to fit the comparison needs
color(1,:) = [0, 0, 0]; % gD
color(2,:) = colorOrig(1,:); % rec 
color(3,:) = colorOrig(5,:); % active
color(4,:) = colorOrig(2,:); % BS
color(5,:) = colorOrig(3,:); % YTE
color(6,:) = colorOrig(9,:); % Gant
color(7,:) = colorOrig(10,:); % Tronti
color(8,:) = colorOrig(1,:); % Rec Chu
colororder(color)

color_plot_blue = ["blue1","blue2","blue3","blue4"];
color_plot_green = ["green1","green2","green3","green4"];
color_plot_red = ["red1","red2","red3","red4","black"];
color_plot_human = ["blue3","green3"];
marker_plot = ["o","+","x","*"];   


%% MAPT

args.groupFit = [250813001:250813004];


%% User input for the runPlot 
time_final = 320; % days
stopTime = hours(days(time_final)); % simulation duration
alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations

%% set up the scenarios
% groups to run
groups = args.groupFit; 

for g = 1:length(groups)
group=groups(g);

% add group parameter if doesn't exist in the model already
p = sbioselect(model,'Name','group','type','parameter');
if isempty(p); p = addparameter(model,'group',0,'Units','dimensionless'); end

% get tasks with the same order as group
[Lia,Locb] = ismember(args.taskTable.Group,group);
tasks = args.taskTable(Lia,:);
[~,p] = sort(Locb(Lia));
tasks = tasks(p,:);

% Display Name from the task spreadsheet
displayName = tasks.DisplayName;

% Create scenarios for the groups
sObj = SimBiology.Scenarios;
add(sObj,'elementwise','group',tasks.Group);

% Add Variants as elementwise to the scenarios
varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
for i= [5,6,7,8] % varCol % Order of columns changed to bring varfit to the right
    add(sObj,'elementwise','variants',tasks.(tasks.Properties.VariableNames{i})(:));
end


% Add doses as elementwise to the scenarios
doseCol = find(contains(lower(tasks.Properties.VariableNames), lower('Dose')));
for i=doseCol
    add(sObj,'elementwise','dose',tasks.(tasks.Properties.VariableNames{i})(:));
end

% Add cartesian variability to the scenarios
if args.variability
    sObj = add_eta_C1s_v0(model,sObj,args);
end

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

    
 
    
         %% Brain protein percent ------------------------------------------------
        

        speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_red(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'CSF'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_red(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'10 mg Q4W','30 mg Q4W','60 mg Q4W','115 mg Q12W'},'Location','southwest','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    
end


%% APP siRNA

args.groupFit = [250904002:250904003];


%% User input for the runPlot 
time_final = 300; % days
stopTime = hours(days(time_final)); % simulation duration
alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations

%% set up the scenarios
% groups to run
groups = args.groupFit; 

for g = 1:length(groups)
group=groups(g);

% add group parameter if doesn't exist in the model already
p = sbioselect(model,'Name','group','type','parameter');
if isempty(p); p = addparameter(model,'group',0,'Units','dimensionless'); end

% get tasks with the same order as group
[Lia,Locb] = ismember(args.taskTable.Group,group);
tasks = args.taskTable(Lia,:);
[~,p] = sort(Locb(Lia));
tasks = tasks(p,:);

% Display Name from the task spreadsheet
displayName = tasks.DisplayName;

% Create scenarios for the groups
sObj = SimBiology.Scenarios;
add(sObj,'elementwise','group',tasks.Group);

% Add Variants as elementwise to the scenarios
varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
for i= [5,6,7,8] % varCol % Order of columns changed to bring varfit to the right
    add(sObj,'elementwise','variants',tasks.(tasks.Properties.VariableNames{i})(:));
end


% Add doses as elementwise to the scenarios
doseCol = find(contains(lower(tasks.Properties.VariableNames), lower('Dose')));
for i=doseCol
    add(sObj,'elementwise','dose',tasks.(tasks.Properties.VariableNames{i})(:));
end

% Add cartesian variability to the scenarios
if args.variability
    sObj = add_eta_C1s_v0(model,sObj,args);
end

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

    
 
    
         %% Brain protein percent ------------------------------------------------
        

        speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,2);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_red(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'CSF'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_red(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'50 mg','75 mg'},'Location','southeast','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    
end



%% HTT -  Group 1

args.groupFit = [250910001:250910002];


%% User input for the runPlot 
time_final = 450; % days
stopTime = hours(days(time_final)); % simulation duration
alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations

%% set up the scenarios
% groups to run
groups = args.groupFit; 

for g = 1:length(groups)
group=groups(g);

% add group parameter if doesn't exist in the model already
p = sbioselect(model,'Name','group','type','parameter');
if isempty(p); p = addparameter(model,'group',0,'Units','dimensionless'); end

% get tasks with the same order as group
[Lia,Locb] = ismember(args.taskTable.Group,group);
tasks = args.taskTable(Lia,:);
[~,p] = sort(Locb(Lia));
tasks = tasks(p,:);

% Display Name from the task spreadsheet
displayName = tasks.DisplayName;

% Create scenarios for the groups
sObj = SimBiology.Scenarios;
add(sObj,'elementwise','group',tasks.Group);

% Add Variants as elementwise to the scenarios
varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
for i= [5,6,7,8] % varCol % Order of columns changed to bring varfit to the right
    add(sObj,'elementwise','variants',tasks.(tasks.Properties.VariableNames{i})(:));
end


% Add doses as elementwise to the scenarios
doseCol = find(contains(lower(tasks.Properties.VariableNames), lower('Dose')));
for i=doseCol
    add(sObj,'elementwise','dose',tasks.(tasks.Properties.VariableNames{i})(:));
end

% Add cartesian variability to the scenarios
if args.variability
    sObj = add_eta_C1s_v0(model,sObj,args);
end

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

    
 
    
         %% Brain protein percent ------------------------------------------------
        

        speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,3);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_red(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'CSF'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_red(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'120 mg Q4W','120 mg Q8W'},'Location','southwest','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:50*24:450*24]);
        xticklabels(ax,{'0','50','100','150','200','250','300','350','400','450'});
        
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    
end


%% HTT -  Group 2

args.groupFit = [250910003:250910007];


%% User input for the runPlot 
time_final = 150; % days
stopTime = hours(days(time_final)); % simulation duration
alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations

%% set up the scenarios
% groups to run
groups = args.groupFit; 

for g = 1:length(groups)
group=groups(g);

% add group parameter if doesn't exist in the model already
p = sbioselect(model,'Name','group','type','parameter');
if isempty(p); p = addparameter(model,'group',0,'Units','dimensionless'); end

% get tasks with the same order as group
[Lia,Locb] = ismember(args.taskTable.Group,group);
tasks = args.taskTable(Lia,:);
[~,p] = sort(Locb(Lia));
tasks = tasks(p,:);

% Display Name from the task spreadsheet
displayName = tasks.DisplayName;

% Create scenarios for the groups
sObj = SimBiology.Scenarios;
add(sObj,'elementwise','group',tasks.Group);

% Add Variants as elementwise to the scenarios
varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
for i= [5,6,7,8] % varCol % Order of columns changed to bring varfit to the right
    add(sObj,'elementwise','variants',tasks.(tasks.Properties.VariableNames{i})(:));
end


% Add doses as elementwise to the scenarios
doseCol = find(contains(lower(tasks.Properties.VariableNames), lower('Dose')));
for i=doseCol
    add(sObj,'elementwise','dose',tasks.(tasks.Properties.VariableNames{i})(:));
end

% Add cartesian variability to the scenarios
if args.variability
    sObj = add_eta_C1s_v0(model,sObj,args);
end

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

    
 
    
         %% Brain protein percent ------------------------------------------------
        

        speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,7);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_red(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'CSF'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_red(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'10 mg Q4W','30 mg Q4W','60 mg Q4W','90 mg Q4W','120 mg Q4W'},'Location','southwest','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    
end

%% Jacifusen

args.groupFit = [250913001,250913003,250913004];


%% User input for the runPlot 
time_final = 450; % days
stopTime = hours(days(time_final)); % simulation duration
alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations

%% set up the scenarios
% groups to run
groups = args.groupFit; 

for g = 1:length(groups)
group=groups(g);

% add group parameter if doesn't exist in the model already
p = sbioselect(model,'Name','group','type','parameter');
if isempty(p); p = addparameter(model,'group',0,'Units','dimensionless'); end

% get tasks with the same order as group
[Lia,Locb] = ismember(args.taskTable.Group,group);
tasks = args.taskTable(Lia,:);
[~,p] = sort(Locb(Lia));
tasks = tasks(p,:);

% Display Name from the task spreadsheet
displayName = tasks.DisplayName;

% Create scenarios for the groups
sObj = SimBiology.Scenarios;
add(sObj,'elementwise','group',tasks.Group);

% Add Variants as elementwise to the scenarios
varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
for i= [5,6,7,8] % varCol % Order of columns changed to bring varfit to the right
    add(sObj,'elementwise','variants',tasks.(tasks.Properties.VariableNames{i})(:));
end


% Add doses as elementwise to the scenarios
doseCol = find(contains(lower(tasks.Properties.VariableNames), lower('Dose')));
for i=doseCol
    add(sObj,'elementwise','dose',tasks.(tasks.Properties.VariableNames{i})(:));
end

% Add cartesian variability to the scenarios
if args.variability
    sObj = add_eta_C1s_v0(model,sObj,args);
end

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

    
 
    
         %% Brain protein percent ------------------------------------------------
        

        speciesModel = 'V_bsf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,4);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_red(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_red(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'Participant 1','Participant 6','Participant 12'},'Location','southwest','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:50*24:450*24]);
        xticklabels(ax,{'0','50','100','150','200','250','300','350','400','450'});
        
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    
end



end
