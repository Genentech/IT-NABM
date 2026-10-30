function runPlotOneSA(model,data,args,options)
% runPlotOneSA  One-way sensitivity sweep of a parameter, plotted as a panel grid.
%   runPlotOneSA(model,data,args, ...) builds scenarios that sweep one parameter
%   (options.p_Name over options.p_Value) at the given dose(s), simulates them,
%   and plots a grid of PK/PD readouts (plasma/brain PK and downstream markers),
%   exporting the figure to args.figLoc.
%
%   Inputs:
%     model   - SimBiology model.
%     data    - dataset table.
%     args    - struct of settings (figure/simulation options).
%     options - name-value options (Group, doses_nmole, p_Name, p_Value,
%               stopTimeDay, SaveName, ...).

arguments
    model SimBiology.Model
    data table
    args struct
    options.Group double = 231000102 % 
    options.DisplayName string = ''
    options.doses_nmole double = [30]*1000/150
    options.p_Name string = 'Km_TfR'
    options.p_Value double = [0.1, 1, 10]
    options.p_DisplayName string = 'TfR Affinity (nM)'
    options.stopTimeDay double = 30
    options.interval_day double = 3000
    options.infusionDur_min double = 0.1
    options.reapeatCount_IV double = 100
    options.SaveName string = 'test'
end

%% User input for the runPlot 
stopTime = hours(days(options.stopTimeDay)); % simulation duration
alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations

%% set up the scenarios

% add group parameter if doesn't exist in the model already
p = sbioselect(model,'Name','group','type','parameter');
if isempty(p); p = addparameter(model,'group',0,'Units','dimensionless'); end

group=options.Group;
% load up the tasksTable and filter for entries for selected groups
tasks = args.taskTable(find(ismember(args.taskTable.Group,group)),:);
% Display Name from the task spreadsheet
if strcmpi(options.DisplayName,'test')
    displayName = tasks.DisplayName;
end

% Create scenarios for the groups
sObj = SimBiology.Scenarios;
add(sObj,'elementwise','group',tasks.Group);
% Add Variants as elementwise to the scenarios
varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
for i=varCol
    add(sObj,'elementwise','variants',tasks.(tasks.Properties.VariableNames{i})(:));
end
LineStyle = [repmat({':'},numel(options.doses_nmole),1);repmat({'-'},numel(options.doses_nmole),1);repmat({'--'},numel(options.doses_nmole),1)];

% Add doses as elementwise to the scenarios
    add(sObj,'elementwise','dose',getdose(model,['ParametericDose_IV']));
    add(sObj,'cartesian','dose_nmole_IV',options.doses_nmole);
    add(sObj,'cartesian',options.p_Name,options.p_Value);


% set up the rest of the parameteric dosing attributes provided to the function
p = sbioselect(model,'Name','stopTimeDay','type','parameter');
p.Value = options.stopTimeDay;
p = sbioselect(model,'Name','interval_day_IV','type','parameter');
p.Value = options.interval_day;
p = sbioselect(model,'Name','infusionDur_min_IV','type','parameter');
p.Value = options.infusionDur_min;
p = sbioselect(model,'Name','reapeatCount_IV','type','parameter');
p.Value = options.reapeatCount_IV;

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
simDataArray = simFunc(sObj, [], [], 0:0.1:stopTime); 

elapsedTime = toc;  % stop the timer and calculate elapsed time
fprintf('The simulation time is %.2f seconds.\n', elapsedTime);

%% Add any post simulation observables that may need the timeseries like min/max
tic
% simDataArray = addobservable(simDataArray, 'Output'     , [char(options.c_Func),'(',char(options.c_Name),'(1:end))'] ,'Units',options.c_Unit);
simDataArray = addobservable(simDataArray, 'dose_mg'    , 'dose_mg_IV(end)'                                          ,'Units','milligram');
        
Obs = vertcat(simDataArray.ScalarObservables);
for i=1:width(Obs)
    scenariosTbl.(Obs.Properties.VariableNames{i}) = Obs.(Obs.Properties.VariableNames{i});
end

elapsedTime = toc;  % stop the timer and calculate elapsed time
fprintf('The observable calculation time is %.2f seconds.\n', elapsedTime);

%% NOTE
% Earlier experiences have shown that addobservable is inefficient for large datasets
% Monitor the time and potentially write your own simple observables

% [ ~, x] = selectbyname(simDataArray, 'V_bsf.AUC_D_nM');
% y= cell2mat(x');
% scenariosTbl.AUC_D_bsf_nMhr = y(end,:)

%% QC
% % Brain PK
% sD = selectbyname(simDataArray, 'V_bsf.D_nM');
% sbioplot(sD)
% % Brain TE
% sD = selectbyname(simDataArray, 'V_bsf.pC1s_tot_norm');
% sbioplot(sD)

%% Plot PK ratios
% create plots 
f = figure; 
f.Position = args.figSize(1,:).*[1,1,1,0.9];
colororder("gem12");
colorOrig = get(gca,'colororder');
color = colorOrig;


% Serum PK
speciesModel = 'V_c.D_tot_ugmL';
    % serum PK
    sp = subplot(3,5,1); hold on
    sp.Box="on";
    
    for i=1:height(scenariosTbl)
    % plot the simulation
    simDataTemp = simDataArray(i);
    simDataTemp = selectbyname(simDataTemp,speciesModel);
    ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{1}),LineStyle=LineStyle{i}, Legend="off");
    end
    % make the plot look nice
    ax.Title.String = options.p_DisplayName;
    ax.FontSize = args.fontSize;
    % y-Axis
    ax.YLim = [1e0, 2e5];
    ax.YScale = 'log';
    ax.YLabel.String = "Plasma PK (ug/mL)";
    % x-Axis
    xTimeAxis(ax,XTimeUnit="day")
    % TASK: write a function to add a second column with different units
    yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);


% Brain PK
speciesModel = 'V_bsf.D_tot_ugmL';
    sp = subplot(3,5,2); hold on
    sp.Box="on";
    
    for i=1:height(scenariosTbl)
    % plot the simulation
    simDataTemp = simDataArray(i);
    simDataTemp = selectbyname(simDataTemp,speciesModel);
    ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{1}),LineStyle=LineStyle{i}, Legend="off");

    end
     % make the plot look nice
    % ax.Title.String = displayName{j};
    ax.FontSize = args.fontSize;
    % y-Axis
    ax.YLim = [0, 4];
    ax.YScale = 'linear'; % linear log
    ax.YLabel.String = "Total Brain PK (ug/mL)";
    % x-Axis
    xTimeAxis(ax,XTimeUnit="day")
    yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);




% C1s free CNS

speciesModel = 'V_bsf.pC1s_ugmL';
    sp = subplot(3,5,3); hold on
    sp.Box="on";
    
    for i=1:height(scenariosTbl)
    % plot the simulation
    simDataTemp = simDataArray(i);
    simDataTemp = selectbyname(simDataTemp,speciesModel);
    ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{1}),LineStyle=LineStyle{i}, Legend="off");

    end
    % make the plot look nice
    % ax.Title.String = displayName{j};
    ax.FontSize = args.fontSize;
    % y-Axis
    ax.YLim = [0, 0.05];
    % ax.YScale = 'log';
    ax.YLabel.String = "Free Brain pC1s (ug/mL)";
    % x-Axis
        xTimeAxis(ax,XTimeUnit="day");
        ax = yAxisSecondUnit(ax,unitLeft='%',unitRight='',MW=80);
    


% C1sComplex Plasma

speciesModel = 'V_bsf.C1sComplex_tot_percent';
    sp = subplot(3,5,4); hold on
    sp.Box="on";
    
    for i=1:height(scenariosTbl)
    % plot the simulation
    simDataTemp = simDataArray(i);
    simDataTemp = selectbyname(simDataTemp,speciesModel);
    ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{1}),LineStyle=LineStyle{i}, Legend="off");

    end
    % make the plot look nice
    % ax.Title.String = displayName{j};
    ax.FontSize = args.fontSize;
    % y-Axis
    ax.YLim = [0, 150];
    % ax.YScale = 'log';
    ax.YLabel.String = "% BL CNS C1sComp";
    % x-Axis
        xTimeAxis(ax,XTimeUnit="day");
        ax = yAxisSecondUnit(ax,unitLeft='%',unitRight='',MW=80);
    



% norm pC4a/C4
speciesModel = 'V_bsf.pC4a_C4_norm';
    sp = subplot(3,5,5); hold on
    sp.Box="on";
    
    for i=1:height(scenariosTbl)
    % plot the simulation
    simDataTemp = simDataArray(i);
    simDataTemp = selectbyname(simDataTemp,speciesModel);
    ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{1}),LineStyle=LineStyle{i}, Legend="off");

    end
    % make the plot look nice
    % ax.Title.String = displayName{j};
    ax.FontSize = args.fontSize;
    % y-Axis
    ax.YLim = [0, 150];
    % ax.YScale = 'log';
    ax.YLabel.String = "Normalized CSF pC4a/C4 (%)";
    % x-Axis
        xTimeAxis(ax,XTimeUnit="day");
        ax = yAxisSecondUnit(ax,unitLeft='%',unitRight='',MW=80);
    


funcName={dbstack().name}; % char(funcName{1})
exportgraphics(f,fullfile(args.figLoc,[char(funcName{1}),'_OneSA_',char(options.SaveName),char(args.figFormat)]),'ContentType','vector','BackgroundColor','none');
