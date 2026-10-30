%% adds variability to input parameters to add bands to the plots

function sObj = add_eta_2TC_v0(model,sObj,args)
% add_eta_2TC_v0  Add parameter-variability samples to a Scenarios object.
%   sObj = add_eta_2TC_v0(model,sObj,args) commits the fitted variant and then
%   adds sampled values for selected variability parameters (e.g. Vm_TfR, eta_CL)
%   to the SimBiology Scenarios object sObj, so that simulations produce
%   variability bands. The sample size is args.varSampleSize.
%
%   Inputs:
%     model - SimBiology model (with the varFit variant).
%     sObj  - SimBiology.Scenarios object to add samples to.
%     args  - struct of settings (varSampleSize).
%   Outputs:
%     sObj  - Scenarios object with variability samples added.

% For random number generator reproducibility
rng(sum(double('jafarnem'))); 

% commit the fit variant ensure the normal distributions are correctly around fitted values and not model defaults
commit(getvariant(model,['varFit']), model);

% add one cartesian element


% p=sbioselect(model,Name='CLp');
% add(sObj,'cartesian','Vplasma',random('Lognormal',log(p.Value),1,[args.varSampleSize,1]))

% For Vplasma as a sample
p=sbioselect(model,Name='Vm_TfR');
add(sObj,'cartesian','Vm_TfR',random('uniform',0.46063,2.237366,[args.varSampleSize,1]))

% Rest if the samples elementwise
p=sbioselect(model,Name='eta_CL');
add(sObj,'elementwise','eta_CL',repelem(random('Lognormal',log(p.Value),.4,[args.varSampleSize,1]),n));



% p=sbioselect(model,Name='eta_1_sig_brain');
% add(sObj,'cartesian','eta_1_sig_brain',random('Lognormal',log(p.Value),.4,[args.varSampleSize,1]));
% 
% % Rest if the samples elementwise
% p=sbioselect(model,Name='eta_1_sig_CSF');
% add(sObj,'elementwise','eta_1_sig_CSF',repelem(random('Lognormal',log(p.Value),.4,[args.varSampleSize,1]),n));


%% Some guidance on using different distributions 
% Available distributions for makedist: https://www.mathworks.com/help/stats/makedist.html

% Most useful distributions are uninorm, normal, loguniform, lognormal
% makedist('Uniform',lowerBound,UpperBound)
% makedist('Normal',mu,sigma) % truncation at zero is needed for most biological parameters
% makedist('Loguniform',lowerBound,UpperBound) % is not avaialble for direct use in random. must create it using make dist
% makedist('Lognormal',mu,sigma)

%% QC framework for distributions

% pd = makedist('uniform',0,100);
% x = random(pd,[10000,1]);
% histogram(x,linspace(-100,200,31)); xscale linear;
% 
% pd = makedist('normal',50,25);
% x = random(pd,[10000,1]);
% histogram(x,linspace(-100,200,31)); xscale linear;
% 
% pd = makedist('normal',50,25);
% tpd = truncate(pd,0,50+3*25);
% x = random(tpd,[10000,1]);
% histogram(x,linspace(-100,200,31)); xscale linear;
% 
% pd = makedist('Loguniform',1,100);
% x = random(pd,[10000,1]);
% histogram(x,logspace(-1,3,21)); xscale log;
% 
% pd = makedist('Lognormal',log(10),1); % sigma of 1 means e-fold (2.7183)
% x = random(pd,[10000,1]);
% histogram(x,logspace(-1,3,21)); xscale log;
% sum((x>10/2.718 & x<10*2.718))/10000*100; % percent within 1 sigma
