function runModel_NOP_Mummery_human(model,data,args)
% runModel_NOP_Mummery_human  Human Mummery clinical scenario for the NOP oligo model.
%   Builds and applies the scenario variants and dosing regimens for the
%   Mummery human study, simulates the model, and plots/saves the PK/PD
%   panels (brain, CSF, plasma oligo, protein and mRNA percent).
%
%   Inputs:
%     model - SimBiology model.
%     data  - dataset table for data overlays.
%     args  - struct of settings.

% Mummery - [250813001:250813004] % 320 120
 % BP41660PET - [250814001]
 % Edwards - [250829001:250829006] % 672
 % Tominersen - [250901001] % 672
 % Test - [250903001:250903002] % Q4W - 300, Q12W - 924, Q6M - 5760
 % BrownKM - Human - [250904001:250904003] % 300
 %% IC50 alone from species
 % BIIB080 Translation mouse and cyno - [250813005:250813006]; 308(Q4W) % 420 (Q12W) % 840 (Q12W - 10 doses)
 % APP siRNA - [250904004]; 1848 (28 days a month); 1800 (30 days a month)
 % Tominersen translation mouse and cyno - [250910011:250910012]; 280
 % Jacifusen - [250913005:250913006]; 308; 300(exactly 28 days after last dose)

 %% All PD from species
 % BIIB080 Translation mouse and cyno - [251301001:251301002]; 308(Q4W) % 420 (Q12W) % 840 (Q12W - 10 doses)
 % APP siRNA - [251401001]; 1848 (28 days a month); 1800 (30 days a month)
 % Tominersen translation mouse and cyno - [251501001:251501002]; 280
 % Jacifusen - [251601001:251601002]; 308; 300(exactly 28 days after last dose)


args.groupFit = [250903001:250903002]; 

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

%% Plot PK
nCol = 2; % 4;
nRow = 2; % 3
nFig = length(group); % ceil(length(group)/nCol);

for j=1:nFig

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

    for i=j % (j-1)*nCol+1:min(j*nCol,length(group))

    %% BSF OligoPK ------------------------------------------------
        
        speciesModel = 'V_b.Oligo_tot_homog_perc'; %'V_c.D_ugmL'; % V_b.Oligo_tot_homog_nM V_bsf.Oligo_tot_homog_nM
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('blue'), DisplayName="BSF Oligo");
       % [~,pk] = peakseek(simDataTemp.Data)
        simDataTemp.Data(end);
        grid_data_sd.mod1(g,1) = trapz(simDataTemp.Time,simDataTemp.Data);
        grid_data_sd.mod2(g,1) = max(simDataTemp.Data);
        
        % Find the relevant data for the scenario
        dataTemp = getData(data, Group=group(i), Analyte={'OligoPK'}, Matrix={'Brain'});
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
        ax.YLim = [0, 100];
        % yticks(ax,[1e2 1e3 1e4 1e5 1e6]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, 300];
        % ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        sgtitle('Oligo in brain (%ID)')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

 %% BDP OligoPK ------------------------------------------------


   
   
        %% CSF OligoPK ------------------------------------------------
        f = figure(2); 
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
        ax.YLim = [1e0, 1e7];
        yticks(ax,[1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        sgtitle('Oligo in CSF (ng/ml)')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on


         %% Plasma OligoPK ------------------------------------------------
        f = figure(3); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_c.Oligo_tot_nM'; %'V_c.D_ugmL';
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
        ax.YLim = [1e-1, 1e7];
        yticks(ax,[1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        sgtitle('Oligo in plasma (ng/ml)')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
         

         %% Brain protein percent ------------------------------------------------
        
 f = figure(4); 
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
        dataTemp = getData(data, Group=group(i), Analyte={'Protein','NoData'}, Matrix={'CSF'});
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
        sgtitle('Protein percent in CSF')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

        
        %% BSF OligoPK ------------------------------------------------
        f = figure(5); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_bsf.Oligo_tot_homog_nM'; %'V_c.D_ugmL'; % V_b.Oligo_tot_homog_nM V_bsf.Oligo_tot_homog_nM
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('blue'), DisplayName="BSF Oligo");
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
        ax.YLim = [1e2, 1e6];
        yticks(ax,[1e2 1e3 1e4 1e5 1e6]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        sgtitle('Oligo in BSF (ng/ml)')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

        %% BDP OligoPK ------------------------------------------------
        
 f = figure(6); 
f.Position = args.figSize(1,:);
        speciesModel = 'V_bdp.Oligo_tot_homog_nM'; %'V_c.D_ugmL'; % V_b.Oligo_tot_homog_nM V_bsf.Oligo_tot_homog_nM
        % serum PK
        % sp = subplot(nRow,nCol,3*floor((g-1)/2)+1);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp = subplot(nRow,nCol,g);     hold on % nRow,nCol,i-(j-1)*nCol); 
        sp.Box="on";
        % plot the simulation
        simDataTemp = simDataArray(find(scenariosTbl.group==group(i)));
        simDataTemp = selectbyname(simDataTemp,speciesModel);
        ax = sbioplotscan(simDataTemp, Name=speciesModel, Parent=sp, Alpha=alpha, Color=colorPicker('blue'), DisplayName="BSF Oligo");
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
        ax.YLim = [1e2, 1e6];
        yticks(ax,[1e2 1e3 1e4 1e5 1e6]);
        % ax.XLim = [stopTime-(168*2 + 24), stopTime];
        ax.XLim = [0, stopTime];
        ax.YScale = 'log';
        ax.YLabel.String = "";
        % x-Axis
        xTimeAxis(ax,XTimeUnit="day")
        sgtitle('Oligo in BDP (nM)')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

        %% Brain protein percent ------------------------------------------------
        
 f = figure(7); 
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
        dataTemp = getData(data, Group=group(i), Analyte={'mRNA'}, Matrix={'BSF'});
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
        sgtitle('mRNA percent in BSF')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on

        %% Brain protein percent ------------------------------------------------
        
 f = figure(8); 
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
        dataTemp = getData(data, Group=group(i), Analyte={'Protein','NoData'}, Matrix={'Cortex'});
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
        sgtitle('Protein percent in BSF')
        % TASK: write a function to add a second column with different units
        % ax = yAxisSecondUnit(ax,unitLeft='ug/mL',unitRight='nM',MW=150);

        hold on
    end

    end
end
end
