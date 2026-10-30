function model = addDose_NOP(model,data,args,options)
% addDose_NOP  Register all study-group and parametric dosing schedules.
%   model = addDose_NOP(model,data,args) adds the named SimBiology dose objects
%   used across the NOP scenarios: a null padding dose, one dose per study
%   group (keyed by a doseGroup<ID> name), and the parametric repeat doses
%   used for automated scenario exploration (IV/SC/IT/ICV routes).
%
%   Inputs:
%     model   - SimBiology model to add doses to.
%     data    - dataset table (study/group metadata).
%     args    - struct of project settings.
%     options - name-value options (Name).
%   Outputs:
%     model   - model with the dose objects added.

arguments
    model SimBiology.Model
    data table
    args struct
    options.Name string = ''
end

%% Null dose that will be used to pad the dose arrays for automated scenario exploration
model = addDose_mg2nmole(model,Name='doseNull'      ,Time_hr=0,Dose_mgkg=0,Infusion_min=0.1,MW=150,Route='IV');


%% Doses for each study group
model = addDose_mg2nmole(model,Name='doseGroup101-100mgCyno-ICV-7kDa-singleBolus'      ,Time_hr=0,Dose_mgkg=100/3.8,Infusion_min=0.1,MW=7  ,SpeciesName='V_cv.Oligo');
model = addDose_mg2nmole(model,Name='doseGroup102-100mgCyno-ICV-14kDa-singleBolus'     ,Time_hr=0,Dose_mgkg=100/3.8,Infusion_min=0.1,MW=14 ,SpeciesName='V_cv.Oligo');
model = addDose_mg2nmole(model,Name='doseGroup103-20mg/kg-IV-102kDa-singleBolus'       ,Time_hr=0,Dose_mgkg=20     ,Infusion_min=0.1,MW=102,SpeciesName='V_c.OligoBrS');
model = addDose_mg2nmole(model,Name='doseGroup104-20mg/kg-IV-110kDa-singleBolus'       ,Time_hr=0,Dose_mgkg=20     ,Infusion_min=0.1,MW=110,SpeciesName='V_c.OligoBrS');

% MAPT ASO Roche EastonA
% Mouse
model = addDose_mg2nmole(model,Name='doseGroup250801001-2mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=2     ,Infusion_min=0.1,MW=6,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250801002-4mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=4     ,Infusion_min=0.1,MW=6,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250801003-6mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=6     ,Infusion_min=0.1,MW=6,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250801004-8mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=8     ,Infusion_min=0.1,MW=6,SpeciesName='V_cv.Oligo_nM');
% Cyno
model = addDose_mg2nmole(model,Name='doseGroup250805001-2x1.18mg/kg-ICV-Cyno'       ,Time_hr=[0,168],Dose_mgkg=1.18     ,Infusion_min=0.1,MW=6,SpeciesName='V_cv.Oligo_nM');

% Biogen ASO
model = addDose_mg2nmole(model,Name='doseGroup250808001-0.12mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=0.12     ,Infusion_min=0.1,MW=6,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250808002-0.36mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=0.36     ,Infusion_min=0.1,MW=6,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250808003-0.6mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=0.6     ,Infusion_min=0.1,MW=6,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250808004-1.2mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=1.2     ,Infusion_min=0.1,MW=6,SpeciesName='V_cv.Oligo_nM');

% Internal BIIB080
model = addDose_mg2nmole(model,Name='doseGroup250811001-0.12mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=0.12     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250811002-0.36mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=0.36     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250811003-1.2mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=1.2     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250811004-3.6mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=3.6     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');

% EastonA - Mouse
model = addDose_mg2nmole(model,Name='doseGroup260602001-100ug-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=0.1/0.025     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');

% Kordasiewicz - Mouse
model = addDose_mg2nmole(model,Name='doseGroup260603001-30ug-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=0.03/0.025     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup260603002-100ug-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=0.1/0.025     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup260603003-500ug-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=0.5/0.025     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');

% Monine
model = addDose_mg2nmole(model,Name='doseGroup250812001-4mg-IT-Cyno'       ,Time_hr=[0,336+1,672+1,1344+1,2016+1],Dose_mgkg=1.05     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250812002-12mg-IT-Cyno'       ,Time_hr=[0,336+1,672+1,1344+1,2016+1],Dose_mgkg=3.16     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250812003-35mg-IT-Cyno'       ,Time_hr=[0,336+1,672+1,1344+1,2016+1],Dose_mgkg=9.21     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250812004-12mg-IT-Cyno'       ,Time_hr=[0,336+1,672+1],Dose_mgkg=3.16     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');


% Monine - copy
model = addDose_mg2nmole(model,Name='doseGroup251208001-4mg-IT-Cyno'       ,Time_hr=[0,336+1,672+1,1344+1,2016+1],Dose_mgkg=1.05     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup251208002-12mg-IT-Cyno'       ,Time_hr=[0,336+1,672+1,1344+1,2016+1],Dose_mgkg=3.16     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup251208003-35mg-IT-Cyno'       ,Time_hr=[0,336+1,672+1,1344+1,2016+1],Dose_mgkg=9.21     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup251208004-12mg-IT-Cyno'       ,Time_hr=[0,336+1,672+1],Dose_mgkg=3.16     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% Mummery
model = addDose_mg2nmole(model,Name='doseGroup250813001-0.14mg/kg-IT-Human'       ,Time_hr=[0,672+1,1344+1,2016+1],Dose_mgkg=0.14     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250813002-0.43mg/kg-IT-Human'       ,Time_hr=[0,672+1,1344+1,2016+1],Dose_mgkg=0.43     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250813003-0.86mg/kg-IT-Human'       ,Time_hr=[0,672+1,1344+1,2016+1],Dose_mgkg=0.86     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250813004-1.64mg/kg-IT-Human'       ,Time_hr=[0,2016+1],Dose_mgkg=1.64     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% BP41660 PET
model = addDose_mg2nmole(model,Name='doseGroup250814001-0.14mg/kg-IT-Human'       ,Time_hr=0,Dose_mgkg=0.14     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% Devos
model = addDose_mg2nmole(model,Name='doseGroup250822001-7.9mg/kg-IT-Cyno'       ,Time_hr=[0],Dose_mgkg=7.9     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250822002-13.16mg/kg-IT-Cyno'       ,Time_hr=[0],Dose_mgkg=13.16     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250822003-2.63-7.9mg/kg-IT-Cyno'       ,Time_hr=[0,168],Dose_mgkg=[2.63,7.9]     ,Infusion_min=[0.1,0.1],MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250822004-2.63-7.9mg/kg-IT-Cyno'       ,Time_hr=[0,168],Dose_mgkg=[2.63,7.9]     ,Infusion_min=[0.1,0.1],MW=7,SpeciesName='V_csf.Oligo_nM');

% Kordasiewicz
model = addDose_mg2nmole(model,Name='doseGroup260529001-4mg-IT-Cyno'       ,Time_hr=[0,14*24,28*24,56*24,84*24],Dose_mgkg=4/3.8     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup260529002-12mg-IT-Cyno'       ,Time_hr=[0,14*24,28*24,56*24,84*24],Dose_mgkg=12/3.8     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% Edwards
model = addDose_mg2nmole(model,Name='doseGroup250829001-Placebo-0.86mg/kg-Q4W-IT-Human'       ,Time_hr=[6048:672:15456],Dose_mgkg=0.86     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250829002-Placebo-1.64mg/kg-Q12W-IT-Human'       ,Time_hr=[6048:2016:14112],Dose_mgkg=1.64     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250829003-0.14-Q4W-0.86mg/kg-Q12W-IT-Human'       ,Time_hr=[0:672:2016,6048:2016:14112],Dose_mgkg=[0.14*ones(1,numel(0:672:2016)),0.86*ones(1,numel(6048:2016:14112))]     ,Infusion_min=0.1*ones(1,numel([0:672:2016,6048:2016:14112])),MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250829004-0.43-Q4W-0.86mg/kg-Q12W-IT-Human'       ,Time_hr=[0:672:2016,6048:2016:14112],Dose_mgkg=[0.43*ones(1,numel(0:672:2016)),0.86*ones(1,numel(6048:2016:14112))]     ,Infusion_min=0.1*ones(1,numel([0:672:2016,6048:2016:14112])),MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250829005-0.86-Q4W-0.86mg/kg-Q12W-IT-Human'       ,Time_hr=[0:672:2016,6048:2016:14112],Dose_mgkg=[0.86*ones(1,numel(0:672:2016)),0.86*ones(1,numel(6048:2016:14112))]     ,Infusion_min=0.1*ones(1,numel([0:672:2016,6048:2016:14112])),MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250829006-1.64-Q12W-1.64mg/kg-Q12W-IT-Human'       ,Time_hr=[0:2016:2016,6048:2016:14112],Dose_mgkg=[1.64*ones(1,numel(0:2016:2016)),1.64*ones(1,numel(6048:2016:14112))]     ,Infusion_min=0.1*ones(1,numel([0:2016:2016,6048:2016:14112])),MW=7,SpeciesName='V_csf.Oligo_nM');

% Torminersen
model = addDose_mg2nmole(model,Name='doseGroup250901001-1.7143mg/kg-IT-Human'       ,Time_hr=[0:672:5376,6048:1344:12768],Dose_mgkg=1.7143     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% BrownKM2021
model = addDose_mg2nmole(model,Name='doseGroup250902001-11.8421mg/kg-IT-Cyno'       ,Time_hr=0,Dose_mgkg=11.8421     ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250902002-15.7895mg/kg-IT-Cyno'       ,Time_hr=0,Dose_mgkg=15.7895     ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250902003-10mg-IT-Cyno'       ,Time_hr=0,Dose_mgkg=10/3.8     ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250902004-15.78956mg/kg-IT-Cyno'       ,Time_hr=0,Dose_mgkg=15.7895     ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250902005-60mg-IT-Cyno'       ,Time_hr=0,Dose_mgkg=60/3.8     ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');

% BrownKM2023
model = addDose_mg2nmole(model,Name='doseGroup260526001-60mg-IT-Cyno'       ,Time_hr=0,Dose_mgkg=60/3.8     ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');

% Test
model = addDose_mg2nmole(model,Name='doseGroup250903001-Xmg/kg-IT-Human'       ,Time_hr=[0:4*7*24:4*7*24*10],Dose_mgkg=115/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250903002-Xmg/kg-IT-Human'       ,Time_hr=[0:4*7*24:4*7*24*10],Dose_mgkg=630/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% BrownKM2021 - Human
model = addDose_mg2nmole(model,Name='doseGroup250904001-0.3571mg/kg-IT-Human'       ,Time_hr=0,Dose_mgkg=0.3571     ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250904002-0.7143mg/kg-IT-Human'       ,Time_hr=0,Dose_mgkg=0.7143     ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250904003-1.0714mg/kg-IT-Human'       ,Time_hr=0,Dose_mgkg=1.0714    ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');

%% HTT
% NorrisDA - Mouse
model = addDose_mg2nmole(model,Name='doseGroup250905001-1.6mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=1.6     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');

% WO_2011_032045_A1_Mouse
model = addDose_mg2nmole(model,Name='doseGroup260525001-50ug-ICV-Mouse'       ,Time_hr=[0:24:336],Dose_mgkg=50e-3/0.025     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');


% TabriziSJ - Cyno
model = addDose_mg2nmole(model,Name='doseGroup250909001-0.5263mg/kg-Q2W-Q4W-IT-Cyno'       ,Time_hr=[0,336,672,1344,2016],Dose_mgkg=0.5263     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250909002-1.8421mg/kg-Q2W-Q4W-IT-Cyno'       ,Time_hr=[0,336,672,1344,2016],Dose_mgkg=1.8421     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250909003-5.2632mg/kg-Q2W-Q4W-IT-Cyno'       ,Time_hr=[0,336,672,1344,2016],Dose_mgkg=5.2632     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250909004-1.8421mg/kg-Q2W-IT-Cyno'       ,Time_hr=[0,336,672],Dose_mgkg=1.8421     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250909005-1.8421mg/kg-IT-Cyno'       ,Time_hr=0,Dose_mgkg=1.8421     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250909006-1.8421mg/kg-Q4W-IT-Cyno'       ,Time_hr=[0,672,1344,2016],Dose_mgkg=1.8421     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250909007-5.2632mg/kg-Q4W-IT-Cyno'       ,Time_hr=[0,672,1344,2016],Dose_mgkg=5.2632     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% Human
model = addDose_mg2nmole(model,Name='doseGroup250910001-1.7143mg/kg-Q4W-IT-Human'       ,Time_hr=[0:672:10080],Dose_mgkg=1.7143     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910002-1.7143mg/kg-Q4W-Q8W-IT-Human'       ,Time_hr=[0,672,2016:1344:10080],Dose_mgkg=1.7143     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910003-0.1429mg/kg-Q4W-IT-Human'       ,Time_hr=[0:672:2016],Dose_mgkg=0.1429     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910004-0.4286mg/kg-Q4W-IT-Human'       ,Time_hr=[0:672:2016],Dose_mgkg=0.4286     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910005-0.8571mg/kg-Q4W-IT-Human'       ,Time_hr=[0:672:2016],Dose_mgkg=0.8571     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910006-1.2857mg/kg-Q4W-IT-Human'       ,Time_hr=[0:672:2016],Dose_mgkg=1.2857     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910007-1.7143mg/kg-Q4W-IT-Human'       ,Time_hr=[0:672:2016],Dose_mgkg=1.7143    ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910008-0.4286mg/kg-Q4W-IT-Human'       ,Time_hr=[0,672],Dose_mgkg=0.4286    ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910009-0.8571mg/kg-Q4W-IT-Human'       ,Time_hr=[0,672],Dose_mgkg=0.8571    ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910010-1.7143mg/kg-Q4W-IT-Human'       ,Time_hr=[0,672],Dose_mgkg=1.7143    ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');


%% Jacifusen
% ShneiderNA - Mouse
model = addDose_mg2nmole(model,Name='doseGroup250911001-12mg/kg-ICV-Mouse'       ,Time_hr=0,Dose_mgkg=12     ,Infusion_min=0.1,MW=7,SpeciesName='V_cv.Oligo_nM');

% ShneiderNA - Cyno
model = addDose_mg2nmole(model,Name='doseGroup250912001-0.7895mg/kg-Q4W-Q8W-IT-Cyno'       ,Time_hr=[0,672,2016,3360,4704,6048],Dose_mgkg=0.7895     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250912002-2.6316mg/kg-Q4W-Q8W-IT-Cyno'       ,Time_hr=[0,672,2016,3360,4704,6048],Dose_mgkg=2.6316     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');


% ShneiderNA - Human
model = addDose_mg2nmole(model,Name='doseGroup250913001-Participant1-IT-Human'       ,Time_hr=[0,294,612,974,1835,2470,3195,3874,4509,5188,5913,6502],Dose_mgkg=[0.2857,0.5714,0.8571,1.1429,1.4286,1.7143,1.7143,1.7143,1.7143,1.7143,1.7143,1.7143]     ,Infusion_min=0.1*ones(1,12),MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250913002-Participant5-IT-Human'       ,Time_hr=[0,340,1019,1654,2334,3014,3670,4328,5710,6525,7046],Dose_mgkg=[1.1429,1.7143,1.7143,1.7143,1.7143,1.7143,1.7143,1.7143,1.7143,1.7143,1.7143]     ,Infusion_min=0.1*ones(1,11),MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250913003-Participant6-IT-Human'       ,Time_hr=[0,634,1314,2016],Dose_mgkg=[1.1429,1.7143,1.7143,1.7143]     ,Infusion_min=0.1*ones(1,4),MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250913004-Participant12-IT-Human'       ,Time_hr=[0,1020,1836,3194,4690],Dose_mgkg=1.7143     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');


%% Translation % Use human weight 70 kg
%% IC50 alone from species
% BIIB080
model = addDose_mg2nmole(model,Name='doseGroup250813005-Xmg/kg-IT-Human'       ,Time_hr=[0:2016:20160],Dose_mgkg=40/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250813006-Xmg/kg-IT-Human'       ,Time_hr=[0:2016:20160],Dose_mgkg=550/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% APP siRNA
model = addDose_mg2nmole(model,Name='doseGroup250904004-Xmg/kg-IT-Human'       ,Time_hr=[0:4320:43200],Dose_mgkg=50/70    ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');

% Tominersen
model = addDose_mg2nmole(model,Name='doseGroup250910011-Xmg/kg-IT-Human'       ,Time_hr=[0:672:6720],Dose_mgkg=200/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250910012-Xmg/kg-IT-Human'       ,Time_hr=[0:672:6720],Dose_mgkg=85/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% Jacifusen
model = addDose_mg2nmole(model,Name='doseGroup250913005-Xmg/kg-IT-Human'       ,Time_hr=[0:672:6720],Dose_mgkg=115/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup250913006-Xmg/kg-IT-Human'       ,Time_hr=[0:672:6720],Dose_mgkg=630/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

%% All PD from species
% BIIB080
model = addDose_mg2nmole(model,Name='doseGroup251301001-Xmg/kg-IT-Human'       ,Time_hr=[0:2016:20160],Dose_mgkg=60/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup251301002-Xmg/kg-IT-Human'       ,Time_hr=[0:2016:20160],Dose_mgkg=700/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% APP siRNA
model = addDose_mg2nmole(model,Name='doseGroup251401001-Xmg/kg-IT-Human'       ,Time_hr=[0:4320:43200],Dose_mgkg=100/70    ,Infusion_min=0.1,MW=13.5,SpeciesName='V_csf.Oligo_nM');

% Tominersen
model = addDose_mg2nmole(model,Name='doseGroup251501001-Xmg/kg-IT-Human'       ,Time_hr=[0:672:6720],Dose_mgkg=450/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup251501002-Xmg/kg-IT-Human'       ,Time_hr=[0:672:6720],Dose_mgkg=60/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% Jacifusen
model = addDose_mg2nmole(model,Name='doseGroup251601001-Xmg/kg-IT-Human'       ,Time_hr=[0:672:6720],Dose_mgkg=115/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup251601002-Xmg/kg-IT-Human'       ,Time_hr=[0:672:6720],Dose_mgkg=630/70     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
%% SA
% model = addDose_mg2nmole(model,Name='doseGroup111111111-10mg-IT-Cyno'       ,Time_hr=[0],Dose_mgkg=2.6316     ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');
model = addDose_mg2nmole(model,Name='doseGroup111111111-115mgQ2W-IT-Human'       ,Time_hr=[0,2016+1],Dose_mgkg=1.64    ,Infusion_min=0.1,MW=7,SpeciesName='V_csf.Oligo_nM');

% % Some Example doses:
% model = addDose_mg2nmole(model,Name='dose-10mg/kg-IV-150kDa-singleBolus'   ,Time_hr=0         ,Dose_mgkg=10          ,Infusion_min=0.1,MW=150,Route='IV');
% model = addDose_mg2nmole(model,Name='dose-10mg/kg-IV-150kDa-QW'            ,Time_hr=24*7*[0:3],Dose_mgkg=10*ones(1,4),Infusion_min=0.1,MW=150,Route='IV');
% model = addDose_mg2nmole(model,Name='dose-10mg/kg-SC-150kDa-singleBolus'   ,Time_hr=0         ,Dose_mgkg=10          ,Infusion_min=0.1,MW=150,Route='SC');

%% Parameteric Oligo Doses

% % Add a parametric mg IV dose
dose  = adddose(model,'DoseParam_mg_Oligo_IV',"repeat");
dose.TargetName            = 'V_c.Oligo_nM';
dose.Amount                = 'dose_Oligo_nmole_IV';  % mg/kg /14*1000 => nmole/kg
dose.StartTime             = 'startTime_day_IV';
dose.Interval              = 'interval_day_IV';
dose.RepeatCount           = 'reapeatCount_IV';
dose.DurationParameterName = 'infusionDur_min_IV';
dose.AmountUnits           = 'nanomole';
dose.TimeUnits             = 'hour';
dose.RateUnits             = 'nanomole/hour';

% % Add a parametric SC dose
dose  = adddose(model,'DoseParam_mg_Oligo_SC',"repeat");
dose.TargetName            = 'V_sc.Oligo_nM';
dose.Amount                = 'dose_Oligo_nmole_SC';
dose.StartTime             = 'startTime_day_SC';
dose.Interval              = 'interval_day_SC';
dose.RepeatCount           = 'reapeatCount_SC';
dose.DurationParameterName = 'infusionDur_min_SC';
dose.AmountUnits           = 'nanomole';
dose.TimeUnits             = 'hour';
dose.RateUnits             = 'nanomole/hour';

% % Add a parametric IT dose
dose  = adddose(model,'DoseParam_mg_Oligo_IT',"repeat");
dose.TargetName            = 'V_csf.Oligo_nM';
dose.Amount                = 'dose_Oligo_nmole_IT';
dose.StartTime             = 'startTime_day_IT';
dose.Interval              = 'interval_day_IT';
dose.RepeatCount           = 'reapeatCount_IT';
dose.DurationParameterName = 'infusionDur_min_IT';
dose.AmountUnits           = 'nanomole';
dose.TimeUnits             = 'hour';
dose.RateUnits             = 'nanomole/hour';


% % Add a parametric ICV dose
dose  = adddose(model,'DoseParam_mg_Oligo_ICV',"repeat");
dose.TargetName            = 'V_cv.Oligo_nM';
dose.Amount                = 'dose_Oligo_nmole_ICV'; 
dose.StartTime             = 'startTime_day_ICV';
dose.Interval              = 'interval_day_ICV';
dose.RepeatCount           = 'reapeatCount_ICV';
dose.DurationParameterName = 'infusionDur_min_ICV';
dose.AmountUnits           = 'nanomole';
dose.TimeUnits             = 'hour';
dose.RateUnits             = 'nanomole/hour';

%% Parameteric OligoBrS Doses

% % Add a parametric mg IV dose
dose  = adddose(model,'DoseParam_mg_OligoBrS_IV',"repeat");
dose.TargetName            = 'V_c.OligoBrS_nM';
dose.Amount                = 'dose_OligoBrS_nmole_IV'; % mg/kg /150*1000 => nmole/kg
dose.StartTime             = 'startTime_day_IV';
dose.Interval              = 'interval_day_IV';
dose.RepeatCount           = 'reapeatCount_IV';
dose.DurationParameterName = 'infusionDur_min_IV';
dose.AmountUnits           = 'nanomole';
dose.TimeUnits             = 'hour';
dose.RateUnits             = 'nanomole/hour';

% % Add a parametric SC dose
dose  = adddose(model,'DoseParam_mg_OligoBrS_SC',"repeat");
dose.TargetName            = 'V_sc.OligoBrS_nM';
dose.Amount                = 'dose_OligoBrS_nmole_SC';
dose.StartTime             = 'startTime_day_SC';
dose.Interval              = 'interval_day_SC';
dose.RepeatCount           = 'reapeatCount_SC';
dose.DurationParameterName = 'infusionDur_min_SC';
dose.AmountUnits           = 'nanomole';
dose.TimeUnits             = 'hour';
dose.RateUnits             = 'nanomole/hour';

% % Add a parametric IT dose
dose  = adddose(model,'DoseParam_mg_OligoBrS_IT',"repeat");
dose.TargetName            = 'V_csf.OligoBrS_nM';
dose.Amount                = 'dose_OligoBrS_nmole_IT';
dose.StartTime             = 'startTime_day_IT';
dose.Interval              = 'interval_day_IT';
dose.RepeatCount           = 'reapeatCount_IT';
dose.DurationParameterName = 'infusionDur_min_IT';
dose.AmountUnits           = 'nanomole';
dose.TimeUnits             = 'hour';
dose.RateUnits             = 'nanomole/hour';


% % Add a parametric ICV dose
dose  = adddose(model,'DoseParam_mg_OligoBrS_ICV',"repeat");
dose.TargetName            = 'V_cv.OligoBrS_nM';
dose.Amount                = 'dose_OligoBrS_nmole_ICV';
dose.StartTime             = 'startTime_day_ICV';
dose.Interval              = 'interval_day_ICV';
dose.RepeatCount           = 'reapeatCount_ICV';
dose.DurationParameterName = 'infusionDur_min_ICV';
dose.AmountUnits           = 'nanomole';
dose.TimeUnits             = 'hour';
dose.RateUnits             = 'nanomole/hour';