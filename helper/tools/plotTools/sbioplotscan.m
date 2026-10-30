function ax = sbioplotscan(simresults, options)
% Creates a percentile plot of scan runs
%
% simresults: array of simdata
%
% Name: names of species to plot as cell array or string array.
% Ex: 'drug' or {'drug','complex'}
%
% Times: vector of double. Timepoints at which al simulation should be
% resampled to be synchronized. Resampling is required to be able to
% compute statistics over all simulations.
% Ex: 0:0.1:100
%
% Parent: axis where to plot
%
% Names: Vector of species or parameter names to plot.
%
% Times: Vector of times to use to compute the statistics.
%
% Alpha: Percent for percentiles between 0 and 1.
%
% Color: Matrix of RGB values between 0 and 1 with one row per species to
% plot.
%
%
% Example: sbioplotscan(simdataObj,Name=["Receptor","Complex"],Alpha=[0.05, 0.01])

arguments
    simresults SimData   {mustBeVector(simresults)}
    options.Name string  {mustBeVector(options.Name),mustBeInSimdata(options.Name,simresults)} = getStateNames(simresults)
    options.Times double {mustBeVector(options.Times)} = getTimes(simresults)
    options.Parent (1,1) {mustBeA(options.Parent,["matlab.ui.control.UIAxes","matlab.graphics.axis.Axes"])}
    options.Alpha (1,:) double {mustBeInRange(options.Alpha,0,1)} = 0.05
    options.Color (:,3) double {mustBeInRange(options.Color,0,1)} = []
    options.YScale (1,1) string {mustBeMember(options.YScale,["linear","log"])} = "linear"
    options.Legend (1,1) string {mustBeMember(options.Legend,["on","off"])} = "on"
    options.DisplayName string  {mustBeVector(options.DisplayName)} = ""
    options.LineStyle string = '-'
    options.LineWidth double = 1.5
end

% Input handling
times = options.Times(:);
speciesname = options.Name(:);
alpha = sort(options.Alpha,'descend');

if ~isempty(options.Color) && height(options.Color)~=numel(speciesname)
    error("Invalid name-value argument 'Color'. " + ...
        "You must specify as many colors as species to plot. " +...
        "Here, the color matrix should contain " + numel(speciesname) + " rows.")
end

% make sure simresults is a column vector
simresults = simresults(:);

% Resample data at tspan
simresults = resample(simresults, times);

% Keep only selected vars
simresults = selectbyname(simresults, speciesname);

% figure setup
if ~isfield(options,'Parent')
    ax = gca;
    cla(ax, 'reset');
else
    ax = options.Parent;
    holdOn = ishold(ax);
end
hold(ax,'on');


% color palette
if isempty(options.Color)

    %     colors = [ 0         0.447        0.741
    %              0.85        0.325        0.098
    %              0.929       0.694        0.125
    %              0.494       0.184        0.556
    %              0.466       0.674        0.188
    %              0.301       0.745        0.933
    %              0.635       0.078        0.184]; % default MATLAB colors

    colors = [  1,0.37,0.41;    % red
        0,0.72,1;       % light blue
        0.39,0.69,0.14; % green
        0.84,0.56,0;    % gold
        0.53,0.5,0.57;  % gray
        0.96,0.40,0.86; % pink
        0.15,0.62,1;    % dark blue
        0,0.74,0.60];   % green blue
else
    colors = options.Color; % user-defined colors
end

oldcolors = colororder(ax);
colororder(ax, colors);

lh = gobjects(numel(speciesname),1);

for i = 1:numel(speciesname)
    % Extract data for specified species
    [ ~, x] = selectbyname(simresults, speciesname(i)) ;
    
    if ~(options.DisplayName=='')
        name = options.DisplayName(i);
    else
        name = speciesname(i);
    end
    
    if strcmpi(class(x),'double')
        x = x;
    else
        x = cell2mat(x');
    end

    % Create Plot
    lh(i) = createPlot(ax,times,x,name,alpha,colors,options);

end


% time units
timeunits = unique(string({simresults(:).TimeUnits}));
if numel(timeunits) == 1
    timeunit = timeunits;
else
    timeunit = "mixed";
end

% data units
datainfos = get(simresults,'DataInfo');
datainfos = vertcat(datainfos{:});
if numel(datainfos) == 1
    dataunits = datainfos.Units;
    dataunit = dataunits;
else
    dataunits = unique(string(cellfun(@(x) x.Units,datainfos,'UniformOutput',false)));
    if numel(dataunits) == 1
        dataunit = dataunits;
    else
        dataunit = "mixed units";
    end
end

% add legend and descriptions
% hold(ax,'off');
% 
if options.Legend == "on"
    lgd = legend;
    lgd.Box = 'off';
end
xlabel(ax,"Time (" + timeunit + ")");
ylabel(ax,dataunit);

% define axes limits
if isMATLABReleaseOlderThan("R2021a")
    offsetFactor = 0.03;
    offset = offsetFactor*times(end);
    ax.XLim = offset*[-1,1]+[0, times(end)];

    minY = min(arrayfun(@(x) min(x.YData),ax.Children));
    maxY = max(arrayfun(@(x) max(x.YData),ax.Children));

    offset = offsetFactor*(maxY-minY);
    ax.YLim = offset*[-1,1]+[minY, maxY];
    %     offset = offsetFactor*range(ax.YLim);
    %     ax.YLim = offset*[-1,1]+ax.YLim;
else
    ax.XLimitMethod = 'padded';
    ax.YLimitMethod = 'padded';
end

if isfield(options,'Parent')

    % restauring hold status if originally 'off'
    if ~holdOn
        hold(ax,'off');
    end

    % restauring original color oder if axes existed
    colororder(ax, oldcolors);

end

ax.YScale = options.YScale;

end % sbioplotscan


%% Code for single Plot
function lh = createPlot(ax,t,y,name,alpha,colors,options)

% sort alpha descending
alpha = sort(alpha, "descend");

% Calculate median
med = median(y, 2);

% line for median first to get color following color order
lh = plot(ax,t, med, 'LineStyle', options.LineStyle, 'LineWidth', options.LineWidth,'Color',colors);
lh.DisplayName = name;
hold(ax,'on');

for i=1:numel(alpha)
    % Calculate percentiles
    pctRange = 100*[alpha(i)/2,1-alpha(i)/2];
    pct = prctile(y, pctRange, 2);

    % percentiles as patch
    xconf = [t; flipud(t)] ;
    yconf = [pct(:,1); flipud(pct(:,2))];

    p = fill(ax, xconf , yconf, 'black');
    p.FaceColor = lh.Color; % same color as line
    p.FaceAlpha = 0.25;
    p.EdgeColor = 'none';
    p.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% change order to get line on top
ax.Children(1:1:(1+numel(alpha))) = ax.Children((1+numel(alpha)):-1:1);

end % createPlot

function times = getTimes(simresults)
% commun data points (should include dosing times)
ref = simresults(1).Time;
for i=2:numel(simresults)
    ref = intersect(ref,simresults(i).time);
end

% interpolate 1000 points between 0 and stoptime
numPoints = 1e3;
timeStart = simresults(1).Time(1);
timeEnd = simresults(1).Time(end);
times = linspace(timeStart, timeEnd, numPoints);

% make sure final vector contains commun time points
times = union(times, ref,'sorted');
end % getTimes

function allnames = getStateNames(simresults)

% get all state names
allnames = string([simresults(1).DataNames]);

% convert the species to compartment.species format to make sure they are
% unique
type = string(cellfun(@(x) x.Type, simresults(1).DataInfo, 'UniformOutput', false));

idxSpecies = type == "species";
if any(idxSpecies)
    speciesInfo  = [simresults(1).DataInfo{idxSpecies}];
    speciesInfo  = struct2table(speciesInfo,"AsArray",true);
    speciesnames = string(speciesInfo.Compartment) + "." + string(speciesInfo.Name);

    % modify species names in list
    allnames(idxSpecies) = speciesnames;
end
end % getStateNames

function mustBeInSimdata(name,simresults)
validnames = getStateNames(simresults);
name = string(name);
tf = ~ismember(name, validnames);

if any(tf)
    eid = 'Names:wrongName';
    msg = "Data for '" + strjoin(name(tf),"', '") + "' not available in SimData objects.";
    throwAsCaller(MException(eid,msg))
end % mustBeInSimdata

end

