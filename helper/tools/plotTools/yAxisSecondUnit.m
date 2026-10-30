function ax = yAxisSecondUnit(ax, options)
% yAxisSecondUnit  Overlay a right-hand y-axis showing a converted unit.
%   ax = yAxisSecondUnit(ax, unitLeft=..., unitRight=..., MW=...) adds a second
%   y-axis on the right whose tick labels show the left-axis values converted to
%   a different unit (e.g. ug/mL <-> nM using the molecular weight MW).
%
%   Inputs:
%     ax      - target axes handle.
%     options - name-value options (unitLeft, unitRight, MW).
%   Outputs:
%     ax      - the (original) left axes handle.

arguments
    ax
    options.unitLeft string = ''
    options.unitRight string = ''
    options.MW double = 150 % must be in kg/mole (Ab: ~150 kg/mole or kDa)
end


ax1 = ax; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'YAxisLocation','right',...
    'Color','none');
% turn off visibility of x axis
ax2.XAxis.Visible = 'off'; 

ax2.YScale = ax1.YScale;
ax2.YLim = ax1.YLim;
ax2.FontSize = ax1.FontSize;

ax1Ytick = ax1.YTick';
ax1.YAxis.TickValues = ax1Ytick;

if ismember(lower(options.unitLeft),lower({'ug/mL','ugmL'}))&&ismember(options.unitRight,{'nM'})
    ax2YTick = ax1Ytick/150*1000;
    ax2.YAxis.TickValues = ax1Ytick;
    ax2.YAxis.TickLabels = compose('%g', round(ax2YTick,2,'significant')); 
    ax2.YAxis.TickLabels{end} = ['(',char(options.unitRight),')']; % BACK
end
if ismember(lower(options.unitLeft),lower({'nM'}))&&ismember(options.unitRight,{'ug/mL','ugmL'})
    ax2YTick = ax1Ytick*150/1000;
    ax2.YAxis.TickValues = ax1Ytick;
    ax2.YAxis.TickLabels = compose('%g', round(ax2YTick,2,'significant')); 
    ax2.YAxis.TickLabels{end} = ['(',char(options.unitRight),')']; % BACK
end
if ismember(lower(options.unitLeft),lower({'%'}))
    ax2YTick = ax1Ytick;
    ax2.YAxis.TickValues = ax1Ytick;
end

    % ax2.YAxis.TickLabels = num2str(ax2YTick,'%g');
    % ax2.YAxis.TickLabels = num2str(round(ax2YTick,2,'significant'),'%g'); % compose('%g', round(ax2YTick,2,'significant')); % BACK

% shrink the x axis a bit to accomodate for the second axis
ax1.Position=ax1.Position.*[1 1 0.8 1];
ax2.Position=ax2.Position.*[1 1 0.8 1];