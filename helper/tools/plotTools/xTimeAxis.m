% automatically take care of the x-axis time label, ranges and ticks
function xTimeAxis(ax, options)
% xTimeAxis  Set a sensible x-axis time label, limits and ticks for a plot.
%   xTimeAxis(ax, XLim=..., XLabel=..., XTick=..., XTimeUnit=...) configures the
%   given axes for time on the x-axis, choosing tick spacing automatically from
%   the axis range and the chosen time unit ('Hour' or 'Day').
%
%   Inputs:
%     ax      - target axes handle.
%     options - name-value options (XLim, XLabel, XTick, XTimeUnit).

arguments
    ax   
    options.XLim double = []
    options.XLabel string = []
    options.XTick double = []
    options.XTimeUnit string = ""
end


% set XLim
if ~isempty(options.XLim)
    ax.XLim = options.XLim;
end


if isempty(options.XTimeUnit)
    options.XTimeUnit = 'Day';
end

switch lower(options.XTimeUnit)
    case 'hour'
        if isempty(options.XTick)
            if ax.XLim(2)>720
                ax.XTick = [0:720:ax.XLim(2)];
                display("It makes more sense to use XTimeUnit = 'Day' in this case");
            elseif (ax.XLim(2)<=720)&(ax.XLim(2)>336) 
                ax.XTick = [0:168:ax.XLim(2)];
            elseif (ax.XLim(2)<=336)&(ax.XLim(2)>168) 
                ax.XTick = [0:48:ax.XLim(2)];
            elseif (ax.XLim(2)<=168)&(ax.XLim(2)>72) 
                ax.XTick = [0:24:ax.XLim(2)];
            elseif (ax.XLim(2)<=72)&(ax.XLim(2)>24) 
                ax.XTick = [0:8:ax.XLim(2)];
            elseif (ax.XLim(2)<=24) 
                ax.XTick = [0:2:ax.XLim(2)];
            end
        else
            ax.XTick = options.XTick;
        end

    case 'day'
        if isempty(options.XTick)
           if ax.XLim(2)>hours(days(75))
                ax.XTick = [0:hours(days(30)):ax.XLim(2)];
            elseif (ax.XLim(2)<=hours(days(75)))&(ax.XLim(2)>hours(days(14))) 
                ax.XTick = [0:hours(days(7)):ax.XLim(2)];
            elseif (ax.XLim(2)<=hours(days(14)))&(ax.XLim(2)>hours(days(7))) 
                ax.XTick = [0:hours(days(2)):ax.XLim(2)];
            elseif (ax.XLim(2)<=hours(days(7))) 
                ax.XTick = [0:hours(days(1)):ax.XLim(2)];
            end
        else
           ax.XTick = hours(days(options.XTick));
        end  
        ax.XTickLabel = num2cell(days(hours(ax.XTick)));
end


% Make axis Tick Length longer for better visibility
ax.TickLength = [0.025, 0.025];

% Changes the scintific notation y-axis into numerical
% ax.YAxis.TickLabels = compose('%g', ax.YAxis.TickValues);

% set Xlabel as Time in Day
ax.XLabel.String = strcat('Time (',lower(options.XTimeUnit),')');
