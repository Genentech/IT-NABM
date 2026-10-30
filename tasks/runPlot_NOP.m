%% Run and plot the H2H molecules single dose

function runPlot_NOP(model,data,args)
% runPlot_NOP  Simulate the fitted groups and plot model-vs-data overlays.
%   runPlot_NOP(model,data,args) loops over the groups in args.groupFit, builds
%   and runs a scenario for each (applying the fitted variant and doses), and
%   plots simulated PK/PD readouts (brain/CSF/plasma/spine oligo PK, mRNA and
%   protein percent) overlaid with the matching observed data, saving the figures.
%
%   Inputs:
%     model - SimBiology model.
%     data  - dataset table for the data overlays.
%     args  - struct of settings (groupFit, figure/simulation options).

%% User input for the runPlot
time_final = 150; % days
stopTime = hours(days(time_final)); % simulation duration
alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations

%% set up the scenarios
% groups to run
groups = args.groupFit; % [210015103,231000201,231000202,231000204]; % [210015103,231000201,231000202]; % [231000102,231000203];  % 231000202; % [240014103]; % 240014103:240014107

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

% Update variant
model = updateVarFit_plot(model, group, args);

% Create scenarios for the groups
sObj = SimBiology.Scenarios;
add(sObj,'elementwise','group',tasks.Group);

% Add Variants as elementwise to the scenarios
varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
for i= [5,6,7,8,4] % varCol % Order of columns changed to bring varfit to the right
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

%% Plot PK
nCol = 2; % 4;
nRow = 2; % 3
nFig = length(group); % ceil(length(group)/nCol);

for j=1:nFig

% create plots 
f = figure(2); 
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

    for i=j % (j-1)*nCol+1:min(j*nCol,length(group))

     %% BSF OligoPK ------------------------------------------------
        
        speciesModel = 'V_bsf.Oligo_tot_homog_nM'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('blue'), DisplayName="BDP Oligo");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'OligoPK'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('blue'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [1e0, 5e5];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on


          %% Deep Brain OligoPK ------------------------------------------------
        f = figure(3); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_bdp.Oligo_tot_homog_nM'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('blue'), DisplayName="BDP Oligo");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'OligoPK'}, Matrix={'Hippo'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('blue'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [1e0, 5e5];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

 %% Brain mRNA percent ------------------------------------------------
        
 f = figure(4); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_bsf.mRNA_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('red'), DisplayName="");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'mRNA'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % meanSD indivDots
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('red'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [0, 110];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
        

        %% Deep Brain mRNA percent ------------------------------------------------
        
 f = figure(5); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_bdp.mRNA_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('red'), DisplayName="");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'mRNA'}, Matrix={'Hippo'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('red'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [0, 110];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
        
         %% Brain protein percent ------------------------------------------------
        
 f = figure(6); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_bsf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('green'), DisplayName="");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'Protein'}, Matrix={'Cortex'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker('green'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [0, 110];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

        %% Deep Brain protein percent ------------------------------------------------
        
 f = figure(7); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_bdp.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('green'), DisplayName="");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'Protein'}, Matrix={'Hippo'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=colorPicker('green'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [0, 110];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

        %% CSF OligoPK ------------------------------------------------
        f = figure(8); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_csf.Oligo_tot_nM'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('blue'), DisplayName="CSF Oligo");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'OligoPK'}, Matrix={'CSF'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('blue'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [5e-2, 1e6];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

         %% Plasma OligoPK ------------------------------------------------
        f = figure(9); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_c.Oligo_tot_ngml'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('blue'), DisplayName="Plasma Oligo");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'OligoPK'}, Matrix={'Plasma'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('blue'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [0, 300];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

         %% Spine OligoPK ------------------------------------------------
        f = figure(10); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_scord.Oligo_tot_homog_ngml'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('blue'), DisplayName="Spine Oligo");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'OligoPK'}, Matrix={'Spine'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('blue'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [1e0, 3e6];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on


        %% CSF protein percent ------------------------------------------------
        
 f = figure(11); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_csf.Prot_perc'; %'V_c.D_ugmL';
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('green'), DisplayName="");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'Protein'}, Matrix={'CSF'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker('green'), LineStyle="none");
        end

        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        pbaspect([2 1 1]);
        grid on; box on;
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        legend off;
        % y-Axis
        ax.YLim = [0, 110];
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    end

    end
end
end
