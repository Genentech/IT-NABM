function ax = plotData(data, options)
% Creates Time vs Conc data plots from the provided and previosuly filter
% data table
% 
% meanSD: calc and plot mean and SD of the timeseries data
% 
% meanError
% 
% meanErrorInput
% 
% geomeanGeoSD
% 
% medianPerc
% 
% indivDots
% 
% percLines
% 
% percBand


% mean/SD, mean/SDInput, geomean/geoSD, median/prctle/%,
% AllIndvidualDots, rangeLines, rangeBand, 

arguments
    data table   
    options.plotType string = "meanSD"
    options.Parent (1,1) {mustBeA(options.Parent,["matlab.ui.control.UIAxes","matlab.graphics.axis.Axes"])}
    options.Color (:,3) double {mustBeInRange(options.Color,0,1)} = []
    options.Error double = []
    options.Percentile (1,2) double = [2.5 97.5]
    options.YScale (1,1) string {mustBeMember(options.YScale,["linear","log"])} = "linear"
    options.Legend (1,1) string {mustBeMember(options.Legend,["on","off"])} = "on"
    options.DisplayName string  {mustBeVector(options.DisplayName)} = ""
    options.Marker string = 'o' % o, +, *, ., x, -, |, square, diamond, ^, v, >, <, pentagram, hexagram, none
    options.MarkerSize double = 5
    options.LineStyle string = '-'
    options.LineWidth double = 1
    % Marker, MarkerSize, LineType, Linewidth, addJitter, 
end

% figure setup
if ~isfield(options,'Parent')
    ax = gca;
    cla(ax, 'reset');
else
    ax = options.Parent;
    holdOn = ishold(ax);
end
hold(ax,'on');

% 
% % color palette to match sbioplotscan
% if isempty(options.Color)
% 
%     %     colors = [ 0         0.447        0.741
%     %              0.85        0.325        0.098
%     %              0.929       0.694        0.125
%     %              0.494       0.184        0.556
%     %              0.466       0.674        0.188
%     %              0.301       0.745        0.933
%     %              0.635       0.078        0.184]; % default MATLAB colors
% 
%     colors = [  1,0.37,0.41;    % red
%         0,0.72,1;       % light blue
%         0.39,0.69,0.14; % green
%         0.84,0.56,0;    % gold
%         0.53,0.5,0.57;  % gray
%         0.96,0.40,0.86; % pink
%         0.15,0.62,1;    % dark blue
%         0,0.74,0.60];   % green blue
% else
%     colors = options.Color; % user-defined colors
% end
% 
% oldcolors = colororder(ax);
% colororder(ax, colors);

% get time and uniqueTime
time = data.Time;
timeUnique = unique(data.Time);

% create different plots
switch options.plotType
    case 'meanSD' 
        dataStat = grpstats(data(:,{'Time','Conc'}),'Time',["mean","std"]);
        % Plot the points or line
        p1 = plot(dataStat.Time, dataStat.mean_Conc); 
        p1.Color=options.Color;
        p1.Marker=options.Marker;
        p1.MarkerSize=options.MarkerSize;
        p1.MarkerFaceColor=options.Color;
        p1.LineStyle=options.LineStyle;
        p1.LineWidth=options.LineWidth;
        if (options.DisplayName==""); p1.Annotation.LegendInformation.IconDisplayStyle = 'off'; end

        % Plot error whiskers
        p2 = errorbar(dataStat.Time, dataStat.mean_Conc, dataStat.std_Conc); 
        p2.Marker="none";
        p2.Color=options.Color;
        p2.LineStyle=options.LineStyle;
        p2.LineWidth=options.LineWidth;
        p2.Annotation.LegendInformation.IconDisplayStyle = 'off';

        if options.Legend == "off"
            p1.Annotation.LegendInformation.IconDisplayStyle = 'on';
            lgd = legend;
            lgd.Box = 'off';
        end

    case 'meanError'
        % Plot the points or line
        p1 = plot(data.Time, data.Conc); 
        p1.Color=options.Color;
        p1.Marker=options.Marker;
        p1.MarkerSize=options.MarkerSize;
        p1.MarkerFaceColor=options.Color;
        p1.LineStyle=options.LineStyle;
        p1.LineWidth=options.LineWidth;
        if (options.DisplayName==""); p1.Annotation.LegendInformation.IconDisplayStyle = 'off'; end

        % Plot error whiskers
        p2 = errorbar(data.Time, data.Conc, options.Error); 
        p2.Marker="none";
        p2.Color=options.Color;
        p2.LineStyle=options.LineStyle;
        p2.LineWidth=options.LineWidth;
        p2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    case 'geomeanGeoSD'
        % dataStat = grpGeoStats(data(:,{'Time','Conc'}),'Time','Conc'); % found an easier way to do the same using the function handles
        dataStat = grpstats(data(:,{'Time','Conc'}),'Time',{@(x)geomean(x,1),@(x)geostd(x,1)});
        dataStat = renamevars(dataStat,["Fun1_Conc","Fun2_Conc"],{['geomean_Conc'],['geostd_Conc']});

        % Plot the points or line
        p1 = plot(dataStat.Time, dataStat.geomean_Conc); 
        p1.Color=options.Color;
        p1.Marker=options.Marker;
        p1.MarkerSize=options.MarkerSize;
        p1.MarkerFaceColor=options.Color;
        p1.LineStyle=options.LineStyle;
        p1.LineWidth=options.LineWidth;
        if (options.DisplayName==""); p1.Annotation.LegendInformation.IconDisplayStyle = 'off'; end

        % Plot error whiskers
        yneg = dataStat.geomean_Conc - dataStat.geomean_Conc./dataStat.geostd_Conc;
        ypos = dataStat.geomean_Conc.*dataStat.geostd_Conc - dataStat.geomean_Conc;
        p2 = errorbar(dataStat.Time, dataStat.geomean_Conc, yneg,ypos); 
        p2.Marker="none";
        p2.Color=options.Color;
        p2.LineStyle=options.LineStyle;
        p2.LineWidth=options.LineWidth;
        p2.Annotation.LegendInformation.IconDisplayStyle = 'off';

    case 'medianPerc'
        dataStat = grpstats(data(:,{'Time','Conc'}),'Time',{@(x)median(x,1),@(x)prctile(x,options.Percentile(1)),@(x)prctile(x,options.Percentile(2))});
        dataStat = renamevars(dataStat,["Fun1_Conc","Fun2_Conc","Fun3_Conc"],{['median_Conc'],['lowerPrctile_Conc'],['upperPrctile_Conc']});

        % Plot the points or line
        p1 = plot(dataStat.Time, dataStat.median_Conc); 
        p1.Color=options.Color;
        p1.Marker=options.Marker;
        p1.MarkerSize=options.MarkerSize;
        p1.MarkerFaceColor=options.Color;
        p1.LineStyle=options.LineStyle;
        p1.LineWidth=options.LineWidth;
        if (options.DisplayName==""); p1.Annotation.LegendInformation.IconDisplayStyle = 'off'; end

        % Plot lower prctile
        p2 = plot(dataStat.Time, dataStat.lowerPrctile_Conc); 
        p2.Marker="none";
        p2.Color=options.Color;
        p2.LineStyle='--';
        p2.LineWidth=options.LineWidth*0.75;
        if (options.DisplayName==""); p2.Annotation.LegendInformation.IconDisplayStyle = 'off'; end

        % Plot upper prctile
        p3 = plot(dataStat.Time, dataStat.upperPrctile_Conc); 
        p3.Marker="none";
        p3.Color=options.Color;
        p3.LineStyle='--';
        p3.LineWidth=options.LineWidth*0.75;
        if (options.DisplayName==""); p3.Annotation.LegendInformation.IconDisplayStyle = 'off'; end

    case 'percLines'
        % mostly copy from 'medianPerc'

    case 'percBand'
        % start with 'medianPerc' and add patch similar to sbioplotscan

    case 'indivDots'
        % Plot the points or line
        IDunique = unique(data.ID);
        for i=1:numel(IDunique)
            dataTemp = data(find(data.ID==IDunique(i)),:);
            p1 = plot(dataTemp.Time, dataTemp.Conc); 
            p1.Color=options.Color;
            p1.Marker=options.Marker;
            p1.MarkerSize=options.MarkerSize;
            p1.MarkerFaceColor='none';
            p1.MarkerEdgeColor=options.Color;
            p1.LineStyle=options.LineStyle;
            p1.LineWidth=options.LineWidth;
            if (options.DisplayName==""); p1.Annotation.LegendInformation.IconDisplayStyle = 'off'; end
        end
    otherwise
        warning('Unexpected plot type. No plot created.')
end % switch plotType
