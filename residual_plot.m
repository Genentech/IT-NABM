%% Residual diagnostic plots
% Plots residuals versus prediction from a completed calibration
% (args.calibration.results, expected in the workspace): the first section
% handles a single CSF-protein response, the second a multi-response fit.

%% Single response
close all
id = (~isnan(args.calibration.results.R));
R_new = args.calibration.results.R(id);
Time_new = args.calibration.results.Data.Time(id)/24;
data_new = args.calibration.results.Data.CSF_Protein_nM(id);
sim_new = data_new+R_new;

plot(sim_new,R_new,'o', LineWidth=2)
xlim([0, 105])
ylim([-20, 20])
hold on
yline(0, '--')

%% Multiple response
close all

R_new = args.calibration.results.R;
Time_new = args.calibration.results.Data.Time/24;
data_new = [args.calibration.results.Data.Hippo_mRNA_nM,args.calibration.results.Data.Hippo_Protein_nM,...
    args.calibration.results.Data.Cortex_mRNA_nM,args.calibration.results.Data.Cortex_Protein_nM,...
    args.calibration.results.Data.CSF_Protein_nM];
sim_new = data_new+R_new;

plot(sim_new,R_new,'ob', LineWidth=2)
xlim([0, 105])
ylim([-15, 15])
hold on
yline(0, '--')
