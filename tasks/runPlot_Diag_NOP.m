%% Run and plot the diagnostics

function runPlot_Diag_NOP(model,data,args,options)
% runPlot_Diag_NOP  Simulate a scenario and plot a full diagnostic PK/PD panel grid.
%   runPlot_Diag_NOP(model,data,args, ...) builds a simulation scenario (modes
%   'group', 'doseVariant' or 'doseParamGroup'), runs it, and plots a 3x5 grid of
%   diagnostic PK panels (plasma, superficial/deep brain, spinal cord, CSF, CV,
%   liver, peripheral) and PD panels (mRNA/protein percent), overlaying any
%   matching data, then exports the figure to args.figLoc.
%
%   Inputs:
%     model   - SimBiology model.
%     data    - dataset table for the data overlays.
%     args    - struct of settings (figure/simulation options).
%     options - name-value options (Type, group, stopTimeDay, Variants, Doses,
%               Route, Format, Dose* parameters, displayName).

arguments
    model SimBiology.Model
    data table
    args struct
    options.Type string = 'group' % 'group' or 'doseVariant' or 'doseParamGroup'
    options.group double = 231000102
    options.stopTimeDay double = 30
    options.Variants
    options.Doses
    options.Route string = 'IV' % 'IV' or 'SC' or IT
    options.Format string = 'Oligo' % 'Oligo' or 'OligoBrS'
    options.DoseAmountMpk = 0 % mg/kg
    options.DoseStartTime = 0 % day
    options.DoseIntervalDay = 0 % day
    options.DoseRepeatCount = 0 % 1 means 2 total doses
    options.DoseInfusionDurMin = 0 % min
    options.displayName string = 'test' % 'custom name'
end


%% User input for the runPlot 
stopTime = hours(days(options.stopTimeDay)); % simulation duration
alpha = 0.05; % alpha 0.05 gives median and 95% interval of the simulations

%% set up the scenarios

% add group parameter if doesn't exist in the model already
p = sbioselect(model,'Name','group','type','parameter');
if isempty(p); p = addparameter(model,'group',0,'Units','dimensionless'); end

m = copyobj(model); model = m; % to make sure the model obj is not mutated

% groups to run
switch lower(options.Type)
    case lower('group')
        group=options.group;
        % load up the tasksTable and filter for entries for selected groups
        tasks = args.taskTable(find(ismember(args.taskTable.Group,group)),:);
        % Display Name from the task spreadsheet
        if strcmpi(options.displayName,'test')
            displayName = tasks.DisplayName;
        else
            displayName = options.displayName;
        end
        % Create scenarios for the groups
        sObj = SimBiology.Scenarios;
        add(sObj,'elementwise','group',tasks.Group);
        % Add Variants as elementwise to the scenarios
        varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
        for i=varCol
            add(sObj,'elementwise','variants',tasks.(tasks.Properties.VariableNames{i})(:));
        end
        % Add doses as elementwise to the scenarios
        doseCol = find(contains(lower(tasks.Properties.VariableNames), lower('Dose')));
        for i=doseCol
            add(sObj,'elementwise','dose',tasks.(tasks.Properties.VariableNames{i})(:));
        end

    case lower('doseVariant')
        group =1;
        Variants = options.Variants;
        Doses = options.Doses;
        displayName = options.displayName;

        % Create scenarios for the groups
        sObj = SimBiology.Scenarios;
        add(sObj,'elementwise','group',group);  
        % Add Variants as elementwise to the scenarios
        for i=1:numel(Variants)
            add(sObj,'elementwise','variants',Variants{i});
        end
        % Add doses as elementwise to the scenarios
        for i=1:numel(Doses)
            add(sObj,'elementwise','dose',Doses{i});
        end

    case lower('doseParamGroup')
        group=options.group;
        % load up the tasksTable and filter for entries for selected groups
        tasks = args.taskTable(find(ismember(args.taskTable.Group,group)),:);
        % Display Name from the task spreadsheet
        if strcmpi(options.displayName,'test')
            displayName = tasks.DisplayName;
        else
            displayName = options.displayName;
        end
        % add the relevant parameteric dose
        
        switch lower(options.Format)
            case lower('Oligo')
                switch lower(options.Route)
                    case lower('IV')
                        Doses = getdose(model,'DoseParam_mg_Oligo_IV');
                        % variant for the dose parameters
                        variantObj = addvariant(model,'varDoseParam');
                        addcontent(variantObj, {'parameter', 'dose_mgkg_IV'        , 'Value', options.DoseAmountMpk});
                        addcontent(variantObj, {'parameter', 'startTime_day_IV'    , 'Value', options.DoseStartTime});
                        addcontent(variantObj, {'parameter', 'interval_day_IV'     , 'Value', options.DoseIntervalDay});
                        addcontent(variantObj, {'parameter', 'reapeatCount_IV'     , 'Value', options.DoseRepeatCount});
                        addcontent(variantObj, {'parameter', 'infusionDur_min_IV'  , 'Value', options.DoseInfusionDurMin});
                    case lower('SC')
                        Doses = getdose(model,'DoseParam_mg_Oligo_SC');
                        % variant for the dose parameters
                        variantObj = addvariant(model,'varDoseParam');
                        addcontent(variantObj, {'parameter', 'dose_mgkg_SC'        , 'Value', options.DoseAmountMpk});
                        addcontent(variantObj, {'parameter', 'startTime_day_SC'    , 'Value', options.DoseStartTime});
                        addcontent(variantObj, {'parameter', 'interval_day_SC'     , 'Value', options.DoseIntervalDay});
                        addcontent(variantObj, {'parameter', 'reapeatCount_SC'     , 'Value', options.DoseRepeatCount});
                        addcontent(variantObj, {'parameter', 'infusionDur_min_SC'  , 'Value', options.DoseInfusionDurMin});
                    case lower('IT')
                        Doses = getdose(model,'DoseParam_mg_Oligo_IT');
                        % variant for the dose parameters
                        variantObj = addvariant(model,'varDoseParam');
                        addcontent(variantObj, {'parameter', 'dose_mgkg_IT'        , 'Value', options.DoseAmountMpk});
                        addcontent(variantObj, {'parameter', 'startTime_day_IT'    , 'Value', options.DoseStartTime});
                        addcontent(variantObj, {'parameter', 'interval_day_IT'     , 'Value', options.DoseIntervalDay});
                        addcontent(variantObj, {'parameter', 'reapeatCount_IT'     , 'Value', options.DoseRepeatCount});
                        addcontent(variantObj, {'parameter', 'infusionDur_min_IT'  , 'Value', options.DoseInfusionDurMin});
                    case lower('ICV')
                        Doses = getdose(model,'DoseParam_mg_Oligo_ICV');
                        % variant for the dose parameters
                        variantObj = addvariant(model,'varDoseParam');
                        addcontent(variantObj, {'parameter', 'dose_mgkg_ICV'        , 'Value', options.DoseAmountMpk});
                        addcontent(variantObj, {'parameter', 'startTime_day_ICV'    , 'Value', options.DoseStartTime});
                        addcontent(variantObj, {'parameter', 'interval_day_ICV'     , 'Value', options.DoseIntervalDay});
                        addcontent(variantObj, {'parameter', 'reapeatCount_ICV'     , 'Value', options.DoseRepeatCount});
                        addcontent(variantObj, {'parameter', 'infusionDur_min_ICV'  , 'Value', options.DoseInfusionDurMin});
                end
            case lower('OligoBrS')
                switch lower(options.Route)
                    case lower('IV')
                        Doses = getdose(model,'DoseParam_mg_OligoBrS_IV');
                        % variant for the dose parameters
                        variantObj = addvariant(model,'varDoseParam');
                        addcontent(variantObj, {'parameter', 'dose_mgkg_IV'        , 'Value', options.DoseAmountMpk});
                        addcontent(variantObj, {'parameter', 'startTime_day_IV'    , 'Value', options.DoseStartTime});
                        addcontent(variantObj, {'parameter', 'interval_day_IV'     , 'Value', options.DoseIntervalDay});
                        addcontent(variantObj, {'parameter', 'reapeatCount_IV'     , 'Value', options.DoseRepeatCount});
                        addcontent(variantObj, {'parameter', 'infusionDur_min_IV'  , 'Value', options.DoseInfusionDurMin});
                    case lower('SC')
                        Doses = getdose(model,'DoseParam_mg_OligoBrS_SC');
                        % variant for the dose parameters
                        variantObj = addvariant(model,'varDoseParam');
                        addcontent(variantObj, {'parameter', 'dose_mgkg_SC'        , 'Value', options.DoseAmountMpk});
                        addcontent(variantObj, {'parameter', 'startTime_day_SC'    , 'Value', options.DoseStartTime});
                        addcontent(variantObj, {'parameter', 'interval_day_SC'     , 'Value', options.DoseIntervalDay});
                        addcontent(variantObj, {'parameter', 'reapeatCount_SC'     , 'Value', options.DoseRepeatCount});
                        addcontent(variantObj, {'parameter', 'infusionDur_min_SC'  , 'Value', options.DoseInfusionDurMin});
                    case lower('IT')
                        Doses = getdose(model,'DoseParam_mg_OligoBrS_IT');
                        % variant for the dose parameters
                        variantObj = addvariant(model,'varDoseParam');
                        addcontent(variantObj, {'parameter', 'dose_mgkg_IT'        , 'Value', options.DoseAmountMpk});
                        addcontent(variantObj, {'parameter', 'startTime_day_IT'    , 'Value', options.DoseStartTime});
                        addcontent(variantObj, {'parameter', 'interval_day_IT'     , 'Value', options.DoseIntervalDay});
                        addcontent(variantObj, {'parameter', 'reapeatCount_IT'     , 'Value', options.DoseRepeatCount});
                        addcontent(variantObj, {'parameter', 'infusionDur_min_IT'  , 'Value', options.DoseInfusionDurMin});
                    case lower('ICV')
                        Doses = getdose(model,'DoseParam_mg_OligoBrS_ICV');
                        % variant for the dose parameters
                        variantObj = addvariant(model,'varDoseParam');
                        addcontent(variantObj, {'parameter', 'dose_mgkg_ICV'        , 'Value', options.DoseAmountMpk});
                        addcontent(variantObj, {'parameter', 'startTime_day_ICV'    , 'Value', options.DoseStartTime});
                        addcontent(variantObj, {'parameter', 'interval_day_ICV'     , 'Value', options.DoseIntervalDay});
                        addcontent(variantObj, {'parameter', 'reapeatCount_ICV'     , 'Value', options.DoseRepeatCount});
                        addcontent(variantObj, {'parameter', 'infusionDur_min_ICV'  , 'Value', options.DoseInfusionDurMin});
                end
        end
        % Create scenarios for the groups
        sObj = SimBiology.Scenarios;
        add(sObj,'elementwise','group',group);  
        % Add Variants as elementwise to the scenarios
        varCol = find(contains(lower(tasks.Properties.VariableNames), lower('Variant')));
        for i=varCol
            add(sObj,'elementwise','variants',tasks.(tasks.Properties.VariableNames{i})(:));
        end
        % Add doses as elementwise to the scenarios 
        for i=1:numel(Doses)
            add(sObj,'elementwise','dose',Doses);
        end
        add(sObj,'elementwise','variants',getvariant(model,'varDoseParam'));
end

% Add cartesian variability to the scenarios
if args.variability
    sObj = add_eta_2TC_v0(model,sObj,args);
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

simDataArray = simFunc(sObj, [], [], 0:1:stopTime);

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
f.Position = args.figSize(1,:);
colororder("gem12");
colorOrig = get(gca,'colororder');
color = colorOrig;

i=1;

    %  Free PK -------------------------------------------------
        % serum PK
        sp = subplot(3,5,1); hold on; grid on;
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_c.BrS_tot_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Total BrS");

        speciesModel = 'V_c.Oligo_tot_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Total Oligo");


        speciesModel = 'V_c.OligoBrS_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle=':', DisplayName="Intact Oligo-BrS");
        
        % speciesModel = 'V_c.Oligo_tot_nM';
        % simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        % simDataTemp = selectbyname(simDataTemp,speciesModel);
        % ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Total Drug");
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'Serum','Plasma'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
    
        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [1e0, 1e5];
        ax.YScale = 'log';
        ax.YLabel.String = "Serum PK (nM)";
        % x-Axis
        ax.XLim = [-0.7, 7.7]*24;
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

    

        % Brain PK
        sp = subplot(3,5,2); hold on; grid on
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_bsf.BrS_tot_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Total BrS in Homogenate");

        speciesModel = 'V_bsf.Oligo_tot_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Total Oligo in Homogenate");

        speciesModel = 'V_bsf.OligoCell_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-.', DisplayName="Oligo in cell in Homogenate");
   
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'CEREBELLUM','FRONTAL_CORTEX','HIPPOCAMPUS'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [1e-1, 1e4];
        ax.YScale = 'log';
        ax.YLabel.String = "Superficial Brain PK (nM)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);


     
        % Brain PK
        sp = subplot(3,5,3); hold on; grid on
        sp.Box="on";
        % plot the simulation
        % plot the simulation
        speciesModel = 'V_bdp.BrS_tot_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Total BrS in Homogenate");

        speciesModel = 'V_bdp.Oligo_tot_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Total Oligo in Homogenate");
      
        speciesModel = 'V_bdp.OligoCell_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-.', DisplayName="Oligo in cell in Homogenate");
     
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'CEREBELLUM','FRONTAL_CORTEX','HIPPOCAMPUS'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [1e-1, 1e4];
        ax.YScale = 'log';
        ax.YLabel.String = "Deep Brain PK (nM)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);


        % 
        % 
        % % Brain PK
        % sp = subplot(3,5,4); hold on; grid on
        % sp.Box="on";
        % % plot the simulation
        %   % plot the simulation
        % speciesModel = 'V_b.BrS_tot_homog_nM';
        % simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        % simDataTemp = selectbyname(simDataTemp,speciesModel);
        % ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Total BrS in Homogenate");
        % 
        % speciesModel = 'V_b.Oligo_tot_homog_nM';
        % simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        % simDataTemp = selectbyname(simDataTemp,speciesModel);
        % ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Total Oligo in Homogenate");
        % 
        % % Find the relevant data for the scenario
        % dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'CEREBELLUM','FRONTAL_CORTEX','HIPPOCAMPUS'});
        % if ~isempty(dataTemp)
        %     % plotData: gets the filtered time vs conc data with option for
        %     % different types of plots
        %     % TASK: make a plot XY version of the plotData function
        %     % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
        %     ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        % end
        %  % make the plot look nice
        % % ax.Title.String = displayName{i};
        % ax.FontSize = args.fontSize;
        % % y-Axis
        % ax.YLim = [1e-1, 1e4];
        % ax.YScale = 'log';
        % ax.YLabel.String = "Whole Brain PK (nM)";
        % % x-Axis
        % xTimeAxis(ax,XTimeUnit="day")
        % % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);



        % Spinal cord PK
        sp = subplot(3,5,4); hold on; grid on
        sp.Box="on";
        % plot the simulation
          % plot the simulation
        speciesModel = 'V_scord.BrS_tot_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Total BrS in Homogenate");

        speciesModel = 'V_scord.Oligo_tot_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Total Oligo in Homogenate");

        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'CEREBELLUM','FRONTAL_CORTEX','HIPPOCAMPUS'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [1e-1, 1e4];
        ax.YScale = 'log';
        ax.YLabel.String = "Spinal Cord PK (nM)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

           
        % CSF PK
        sp = subplot(3,5,12); hold on; grid on
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_csf.BrS_tot_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Total BrS");

        speciesModel = 'V_csf.Oligo_tot_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Total Oligo");

        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix='CSF');
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [1e0, 1e5];
        ax.YScale = 'log';
        ax.YLabel.String = "CSF PK (nM)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);


    
        % CV PK
        sp = subplot(3,5,11); hold on; grid on
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_cv.OligoBrS_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Oligo-BrS");

        speciesModel = 'V_cv.Oligo_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Oligo");

        speciesModel = 'V_cv.BrS_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle=':', DisplayName="BrS");

        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix='ICV');
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=sp, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [1e0, 1e5];
        ax.YScale = 'log';
        ax.YLabel.String = "Cerebroventricular PK (nM)";
        % x-Axis
        ax.XLim = [-0.7, 7.7]*24;
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);


     
        % Liver PK
        sp = subplot(3,5,5); hold on; grid on
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_liver.BrS_tot_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Total BrS in Homogenate");

        speciesModel = 'V_liver.Oligo_tot_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Total Oligo in Homogenate");
        
        speciesModel = 'V_liver.OligoCell_homog_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-.', DisplayName="Oligo in cell in Homogenate");
     
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'CEREBELLUM','FRONTAL_CORTEX','HIPPOCAMPUS'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [1e0, 1e5];
        ax.YScale = 'log';
        ax.YLabel.String = "Liver PK (nM)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);


    
        %  Peripheral PK
        sp = subplot(3,5,6); hold on; grid on
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_p.BrS_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Total BrS");

        speciesModel = 'V_p.OligoBrS_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle=':', DisplayName="Intact Oligo-BrS");
 
        speciesModel = 'V_p_Oligo.Oligo_nM';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Free Oligo (separate comp)");
               
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'Serum','Plasma'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
    
        % make the plot look nice
        ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [1e0, 1e5];
        ax.YScale = 'log';
        ax.YLabel.String = "Peripheral PK (nM)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % TASK: write a function to add a second column with different units
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);


        %--------------PD-----------------------

        % Brain PD
        sp = subplot(3,5,7); hold on; grid on
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_bsf.mRNA_perc';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="mRNA");
     
        speciesModel = 'V_bsf.Prot_perc';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Protein");

        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'CEREBELLUM','FRONTAL_CORTEX','HIPPOCAMPUS'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [0, 120];
        ax.YScale = 'linear';
        ax.YLabel.String = "Superficial Brain PD (%)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);


        % Brain PD
        sp = subplot(3,5,8); hold on; grid on
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_bdp.mRNA_perc';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="mRNA");
     
        speciesModel = 'V_bdp.Prot_perc';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Protein");

        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'CEREBELLUM','FRONTAL_CORTEX','HIPPOCAMPUS'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [0, 120];
        ax.YScale = 'linear';
        ax.YLabel.String = "Deep Brain PD (%)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);



        % Spinal cord PD
        sp = subplot(3,5,9); hold on; grid on
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_scord.mRNA_perc';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="mRNA");

        speciesModel = 'V_scord.Prot_perc';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Protein");

        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'scord'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [0, 120];
        ax.YScale = 'linear';
        ax.YLabel.String = "Spinal cord PD (%)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);



        % Liver PD
        sp = subplot(3,5,10); hold on; grid on
        sp.Box="on";
        % plot the simulation
        speciesModel = 'V_liver.mRNA_perc';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="mRNA");

        speciesModel = 'V_liver.Prot_perc';
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=ax, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='--', DisplayName="Protein");

        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'liver'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [0, 120];
        ax.YScale = 'linear';
        ax.YLabel.String = "Liver PD (%)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);


        % CSF PD
        sp = subplot(3,5,13); hold on; grid on
        sp.Box="on";
        % plot the simulation
        % speciesModel = 'V_csf.Prot_perc';
        % simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        % simDataTemp = selectbyname(simDataTemp,speciesModel);
        % ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker(tasks.Color{i}), LineStyle='-', DisplayName="Protein");

        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte='PK', Matrix={'CEREBELLUM','FRONTAL_CORTEX','HIPPOCAMPUS'});
        if ~isempty(dataTemp)
            % plotData: gets the filtered time vs conc data with option for
            % different types of plots
            % TASK: make a plot XY version of the plotData function
            % axd = plotData(dataTemp, plotType='meanSD', Parent=sp, Color=color(j,:), LineStyle="none");
            ax = plotData(dataTemp, plotType='indivDots', Parent=ax, Color=colorPicker(tasks.Color{i}), LineStyle="none");
        end
         % make the plot look nice
        % ax.Title.String = displayName{i};
        ax.FontSize = args.fontSize;
        % y-Axis
        ax.YLim = [0, 120];
        ax.YScale = 'linear';
        ax.YLabel.String = "CSF Protein (%)";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        % yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);




%% saving the whole figure here
% for some builds it makes sense to save individual axis separately and
% also build on them. exportgraphics can take ax as input and the naming
% can be adjustes to include specific treatment, dose, matrix or analyte
% infomration of interest
% if want to add date: [char(datetime("today"),'yyyyMMdd_'),dbstack().name,'_']
% funcName = char(dbstack().name);
% append(dbstack().name,'_test',args.figFormat)
funcName={dbstack().name}; % char(funcName{1})
exportgraphics(f,fullfile(args.figLoc,append(char(funcName{1}),'_',num2str(group),'_diagnostics',args.figFormat)),'ContentType','vector','BackgroundColor','none');
