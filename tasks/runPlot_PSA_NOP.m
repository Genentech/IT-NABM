function runPlot_PSA_NOP(args)
% Plot PRCC and/or Sobol global sensitivity results and save the figures.
% Standalone use: call runPlot_PSA_NOP() with no args to load the saved
% results from calibration/PSA_results_NOP.mat (no analysis re-run needed).

arguments
    args struct = loadPSAResults()
end

% ---- PRCC: grouped horizontal bar of rho; '*' marks p >= 0.05 ----
if isfield(args,'PRCC') && isfield(args.PRCC,'rho')
    rho        = args.PRCC.rho;
    pval       = args.PRCC.pval;
    paramNames = args.PRCC.inputNamesDisplay;
    outNames   = args.PRCC.outputNamesDisplay;

    f = figure;
    f.Position = args.figSize(1,:);
    h = barh(rho);
    ax = gca;
    ax.YTick      = 1:numel(paramNames);
    ax.YTickLabel = paramNames;
    ax.XLim       = [-1 1];
    ax.XLabel.String = 'PRCC';
    ax.FontSize   = args.fontSize;
    title('PRCC Global Sensitivity');
    legend(outNames,'Location','bestoutside');
    grid on; box on;

    for i = 1:numel(outNames)
        for j = 1:numel(paramNames)
            if pval(j,i) >= 0.05
                text(h(i).XEndPoints(j), h(i).YEndPoints(j), ' *', ...
                    'FontSize', args.fontSize, 'VerticalAlignment','middle');
            end
        end
    end

    exportgraphics(f, fullfile(args.figLoc,['runPlot_PSA_NOP_PRCC',char(args.figFormat)]), ...
        'ContentType','vector','BackgroundColor','none');
end

% ---- Sobol: first-order + interaction via the SobolResults bar method ----
% bar(SobolResults) creates its OWN figure and returns that figure handle
% (see script_NOP_calibration.m), so size/export the returned handle.
if isfield(args,'Sobol') && isfield(args.Sobol,'results')
    f = bar(args.Sobol.results);
    f.Position = args.figSize(1,:);
    set(findall(f,'-property','FontSize'),'FontSize',args.fontSize);

    exportgraphics(f, fullfile(args.figLoc,['runPlot_PSA_NOP_Sobol',char(args.figFormat)]), ...
        'ContentType','vector','BackgroundColor','none');
end

end

function args = loadPSAResults()
% Load saved PSA results. The figure settings (figSize/figFormat/figLoc) are
% re-read live from args_NOP.xlsx so they are not frozen in the .mat snapshot.
here     = fileparts(mfilename('fullpath'));                   % /tasks
projRoot = fileparts(here);                                    % project root
S = load(fullfile(projRoot,'calibration','PSA_results_NOP.mat'),'PSA');
args = S.PSA;

% Refresh figure output settings from args_NOP.xlsx
argsLoc = fullfile(projRoot,'args_NOP.xlsx');
if isfile(argsLoc)
    a = loadArgs(argsLoc);
    args.figSize   = a.figSize;
    args.figFormat = a.figFormat;
    args.figLoc    = fullfile(projRoot, a.figLoc);  % absolute (xlsx value is relative)
    if ~isfolder(args.figLoc); mkdir(args.figLoc); end
end
end
