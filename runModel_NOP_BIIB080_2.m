function runModel_NOP_BIIB080_2(model,data,args)
% runModel_NOP_BIIB080_2  Run and plot the BIIB080 (tau ASO) cross-species scenario.
%   Builds and applies the variants and doses for the BIIB080 scenario across
%   mouse, cyno, and human groups, simulates the model for each group, and
%   plots/saves the brain/CSF PK and PD (Oligo, mRNA, protein) panels.
%
%   Inputs:
%     model - SimBiology model.
%     data  - dataset table for data overlays.
%     args  - struct of settings.
%% Plot PK
nCol = 4; % 4;
nRow = 4; % 3
nFig = 9; % ceil(length(group)/nCol);



% create plots 
f = figure; 
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


%% MOUSE PK

args.groupFit = [260602001];


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
for i= [5,6,7,8,9] % varCol % Order of columns changed to bring varfit to the right
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

    %% Brain OligoPK ------------------------------------------------
        
        speciesModel = 'V_bsf.Oligo_tot_homog_nM'; %'V_c.D_ugmL'; % V_b.Oligo_tot_homog_nM V_bsf.Oligo_tot_homog_nM
        % serum PK
        sp = subplot(nRow,nCol,1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_blue(g)), DisplayName="BSF Oligo");
       
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'OligoPK'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_blue(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        % legend on;
        legend({'100 ug'},'Location','northeast','FontSize',10,'FontWeight','bold','IconColumnWidth',15);
        title('Mouse');
        % y-Axis
        ax.YLim = [1e0, 5e4];
        yticks(ax,[1e0 1e1 1e2 1e3 1e4]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ylabel({'Brain PK (nM)','(Cortex)'});
        ax.YLabel.FontWeight = 'bold';
        ax.YLabel.FontSize = 14;
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        ax.XLabel.String = "";
        % sgtitle('Oligo in brain (nM)')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

 
   
    
end


%% MOUSE PD

args.groupFit = [260603001:260603003];


%% User input for the runPlot 
time_final = 120; % days
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
for i= [5,6,7,8,9] % varCol % Order of columns changed to bring varfit to the right
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

    

    %% Brain mRNA percent ------------------------------------------------
        

        speciesModel = 'V_bsf.mRNA_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,5);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_blue(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'mRNA'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_blue(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'30 ug','100 ug','300 ug'},'Location','northwest','FontSize',10,'FontWeight','bold','IconColumnWidth',15);
        

        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ylabel({'Brain mRNA','(% Baseline)','(Cortex)'});
        ax.YLabel.FontWeight = 'bold';
        ax.YLabel.FontSize = 14;

        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

        %% Brain protein percent ------------------------------------------------
        

        speciesModel = 'V_bsf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,9);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_blue(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_blue(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        % legend({'0.5 nmol','1.5 nmol','5.0 nmol','15.0 nmol'},'Location','southeast','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ylabel({'Brain Protein','(% Baseline)','(Cortex)'});
        ax.YLabel.FontWeight = 'bold';
        ax.YLabel.FontSize = 14;
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    

         %% CSF protein percent ------------------------------------------------
        

        speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,13);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_blue(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'CSF'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_blue(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        % legend({'0.5 nmol','1.5 nmol','5.0 nmol','15.0 nmol'},'Location','southeast','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off

        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ylabel({'CSF Protein','(% Baseline)'});
        ax.YLabel.FontWeight = 'bold';
        ax.YLabel.FontSize = 14;
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        ax.XLabel.String = "Time (day)";
        ax.XLabel.FontSize = 14;
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    
end



%% CYNO

args.groupFit = [250812001:250812004];


%% User input for the runPlot 
time_final = 180; % days
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
for i= [5,6,7,8,9] % varCol % Order of columns changed to bring varfit to the right
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

    %% BSF OligoPK ------------------------------------------------
        
        speciesModel = 'V_bsf.Oligo_tot_homog_nM'; %'V_c.D_ugmL'; % V_b.Oligo_tot_homog_nM V_bsf.Oligo_tot_homog_nM
        % serum PK
        sp = subplot(nRow,nCol,2);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_green(g)), DisplayName="BSF Oligo");
       
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'OligoPK'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_green(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'4 mg Q2W-Q4W','12 mg Q2W-Q4W','35 mg Q2W-Q4W','12 mg Q2W'}, ...
            'Location','south','NumColumns',2,'FontSize',10, ...
            'FontWeight','bold','IconColumnWidth',15);
        title('Cyno');

        % y-Axis
        ax.YLim = [1e0, 5e4];
        yticks(ax,[1e0 1e1 1e2 1e3 1e4]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        ax.XLabel.String = "";
        % sgtitle('Oligo in brain (nM)')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

       
    
end

%% Cyno PD

args.groupFit = [260529001:260529002];


%% User input for the runPlot 
time_final = 180; % days
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
for i= [5,6,7,8,9] % varCol % Order of columns changed to bring varfit to the right
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


%% Brain mRNA percent ------------------------------------------------
        

        speciesModel = 'V_bsf.mRNA_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,6);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_green(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'mRNA'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_green(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'4mg Q2W-Q4W','12mg Q2W-Q4W'},'Location','southeast','FontSize',10,'FontWeight','bold','IconColumnWidth',15);
        
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
    


%% Brain protein percent ------------------------------------------------
        

        speciesModel = 'V_bsf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,10);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_green(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_green(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        % legend({'10 mg - 30 mg'},'Location','southeast','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off

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
    


        %% CSF protein percent ------------------------------------------------
        

        speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,14);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_green(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'CSF'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_green(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        % legend({'10 mg - 30 mg'},'Location','southeast','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        ax.XLabel.String = "Time (day)";
        ax.XLabel.FontSize = 14;
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    

       
       
    
end


%% HUMAN TRANSLATION

args.groupFit = [251301001:251301002];


%% User input for the runPlot 
time_final = 840; % days
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
for i= [5,6,7,8,9] % varCol % Order of columns changed to bring varfit to the right
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

    %% BSF OligoPK ------------------------------------------------
        
        speciesModel = 'V_bsf.Oligo_tot_homog_nM'; %'V_c.D_ugmL'; % V_b.Oligo_tot_homog_nM V_bsf.Oligo_tot_homog_nM
        % serum PK
        sp = subplot(nRow,nCol,4);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_human(g)), DisplayName="BSF Oligo");
       
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'OligoPK'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_human(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'60 mg Q12W','700 mg Q12W'},'Location','southeast','FontSize',10,'FontWeight','bold','IconColumnWidth',15);
        title('Human prediction');

        % y-Axis
        ax.YLim = [1e0, 1e5];
        yticks(ax,[1e0 1e1 1e2 1e3 1e4]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % xticks(ax,[0:2400:16800])
        ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:168*24:840*24])
        xticklabels(ax,{'0','168','336','504','672','840'}); 
        ax.XLabel.String = "";
        % sgtitle('Oligo in brain (nM)')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

        %% Brain mRNA percent ------------------------------------------------
        

        speciesModel = 'V_bsf.mRNA_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,8);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_human(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_human(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        % legend({'4 mg Q4W','150 mg Q4W'},'Location','southeast','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:168*24:840*24])
        xticklabels(ax,{'0','168','336','504','672','840'}); 
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on


        %% Brain protein percent ------------------------------------------------
        

        speciesModel = 'V_bsf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,12);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_human(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_human(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        % legend({'4 mg Q4W','150 mg Q4W'},'Location','southeast','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:168*24:840*24])
        xticklabels(ax,{'0','168','336','504','672','840'}); 
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on


         %% CSF protein percent ------------------------------------------------
        

        speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,16);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_human(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'CSF'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_human(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        % legend({'4 mg Q4W','150 mg Q4W'},'Location','southeast','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:168*24:840*24])
        xticklabels(ax,{'0','168','336','504','672','840'}); 
        ax.XLabel.String = "Time (day)";
        ax.XLabel.FontSize = 14;
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
end

%% HUMAN DATA
% args.groupFit = [250813001:250813004];
% 
% %% User input for the runPlot 
% time_final = 320; % days
% stopTime = hours(days(time_final)); % simulation duration
% alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations
% 
% %% set up the scenarios
% % groups to run
% groups = args.groupFit; 
% 
% for g = 1:length(groups)
% group=groups(g);
% 
% 
%          %% Brain protein percent ------------------------------------------------
% 
% 
%         speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
%         % serum PK
%         % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
%         sp = subplot(nRow,nCol,12);     hold on % nRow,nCol,i-(j-1)*nCol); 
%         sp.Box="on";
% 
% 
%         % Find the relevant data for the scenario
%         dataTemp = getData(data, Group=groups(g), Analyte={'Protein'}, Matrix={'CSF'});
%         if ~isempty(dataTemp)
%             % plotData: gets the filtered time vs conc data with option for
%             % different types of plots
%             % TASK: make a plot XY version of the plotData function
%             % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
%             ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_red(g)), LineStyle="-",Legend="off");
%         end
% 
%         % make the plot look nice
%         % ax.Title.String = displayName{i};
%         ax.FontSize = args.fontSize;
%         pbaspect([2 1 1]);
%         grid on; box on;
%         set(gca,'XMinorGrid','on');
%         set(gca,'YMinorGrid','on');
%         legend({'10 mg Q4W','30 mg Q4W','60 mg Q4W','115 Q12W'},'Location','southwest','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
% 
%         % y-Axis
%         ax.YLim = [0, 110];
%         yticks(ax,[0 25 50 75 100]);
%         % ax.XLim = [stopTime-(168*2 + 24), stopTime];
%         ax.XLim = [0, stopTime];
%         % ax.YScale = 'log';
%         ax.YLabel.String = "";
%         ax.XLabel.String = "";
%         % x-Axis
%         xTimeAxis(ax,XTimeUnit="day")
%         ax.XLabel.String = "";
%         % sgtitle('Tau protein percent in brain')
%         % TASK: write a function to add a second column with different units
%         % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);
% 
%         hold on
% end



%% Human PD actual

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
for i= [5,6,7,8,9] % varCol % Order of columns changed to bring varfit to the right
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

    
 
    
         %% CSF protein percent ------------------------------------------------
        

        speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,15);     hold on % nRow,nCol,i-(j-1)*nCol); 
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
        % legend({'10 mg Q4W','30 mg Q4W','60 mg Q4W','115 mg Q12W'},'Location','southwest','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:1440:7200])
        xticklabels(ax,{'0','60','120','180','240','300'}); 
        ax.XLabel.String = "Time (day)";
        ax.XLabel.FontSize = 14;
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    


        %% Brain protein percent ------------------------------------------------
        

        speciesModel = 'V_bsf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,11);     hold on % nRow,nCol,i-(j-1)*nCol); 
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
        % legend({'10 mg Q4W','30 mg Q4W','60 mg Q4W','115 mg Q12W'},'Location','southwest','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:1440:7200])
        xticklabels(ax,{'0','60','120','180','240','300'}); 
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

        %% Brain mRNA percent ------------------------------------------------
        

        speciesModel = 'V_bsf.mRNA_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,7);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_red(g)), DisplayName="");
     
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'mRNA'}, Matrix={'Cortex'});
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
        % legend({'10 mg Q4W','30 mg Q4W','60 mg Q4W','115 mg Q12W'},'Location','southwest','FontSize',6,'FontWeight','bold','IconColumnWidth',15);
        legend off
        % y-Axis
        ax.YLim = [0, 110];
        yticks(ax,[0 25 50 75 100]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:1440:7200])
        xticklabels(ax,{'0','60','120','180','240','300'}); 
        ax.XLabel.String = "";
        % sgtitle('Tau protein percent in brain')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on


        %% BSF OligoPK ------------------------------------------------
        
        speciesModel = 'V_bsf.Oligo_tot_homog_nM'; %'V_c.D_ugmL'; % V_b.Oligo_tot_homog_nM V_bsf.Oligo_tot_homog_nM
        % serum PK
        sp = subplot(nRow,nCol,3);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==groups(g)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(color_plot_red(g)), DisplayName="BSF Oligo");
       
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=groups(g), Analyte={'OligoPK'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker(color_plot_human(g)), LineStyle="none");
        end

        % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid off; box on;
        % set(gca,'XMinorGrid','on');
        % set(gca,'YMinorGrid','on');
        legend({'10 mg Q4W','30 mg Q4W','60 mg Q4W','115 mg Q12W'}, ...
            'Location','northeast','NumColumns',2,'FontSize',10, ...
            'FontWeight','bold','IconColumnWidth',15);
        title('Human')
        % y-Axis
        ax.YLim = [1e0, 5e4];
        yticks(ax,[1e0 1e1 1e2 1e3 1e4]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        xticks(ax,[0:1440:7200])
        ax.YScale = 'log';
        ax.YLabel.String = "";
        
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        xticks(ax,[0:1440:7200])
        xticklabels(ax,{'0','60','120','180','240','300'}); 
        ax.XLabel.String = "";
        % sgtitle('Oligo in brain (nM)')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
end

append(dbstack().name,'',args.figFormat)
funcName={dbstack().name};
exportgraphics(f,fullfile(args.figLoc,append(char(funcName{1}),args.figFormat)),'ContentType','vector','BackgroundColor','none');

end
