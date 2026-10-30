%% Build a model tool model to test avidity

function model = buildModel_NOP(args)
% buildModel_NOP  Construct the Neuro Oligo Platform (NOP) SimBiology model.
%   model = buildModel_NOP(args) builds and returns a SimBiology model object
%   describing oligonucleotide PK/PD across the plasma, peripheral, brain
%   (superficial/deep), CSF, spinal cord and liver compartments. It adds all
%   compartments, species, observables, parameters, reactions and the
%   parametric dose definitions, then returns the assembled model.
%
%   Inputs:
%     args  - struct of project settings/paths (see loadArgs).
%   Outputs:
%     model - configured SimBiology model object.

% TASKs:
% 1. Add a thrid brain region of interest for substantia nigra or any
% other future disease-spefcific regions of interest
% 2. Update the diffusion to be molecule size dependent 
% 3. Add Tau protein basline dynamics




% Model Object
model = sbiomodel('NeuroOligoPlatform');
config = getconfigset(model);
options = get(config,'CompileOptions');
set(options,'UnitConversion',true);
set(config,'TimeUnits','hour');
set(config,'SolverType','ode15s'); % Model Settings ode15s
start_time = 0.0; % solution start time [days]
time_step = 0.1; % solution output interval [days] (0.01 days ~ 15 mins)
end_time = 15 *3; % solution end time [days]
time = start_time:time_step:end_time;
set(config.SolverOptions,'OutputTimes',time);
set(config,'StopTime',time(end));
tol_abs = 1e-6;
tol_rel = 1e-3;
set(config.SolverOptions,'AbsoluteTolerance',tol_abs);
set(config.SolverOptions,'RelativeTolerance',tol_rel);

% Define accessory units  
addparameter(model,'day' ,1,'ValueUnits','day','ConstantValue',true);
addparameter(model,'kg'  ,1,'ValueUnits','kilogram','ConstantValue',true);
addparameter(model,'nM'  ,1,'ValueUnits','nanomolarity','ConstantValue',true);
addparameter(model,'mL'  ,1,'ValueUnits','milliliter','ConstantValue',true);
addparameter(model,'ngml',1,'ValueUnits','nanogram/milliliter','ConstantValue',true);
addparameter(model,'ugmL',1,'ValueUnits','microgram/milliliter','ConstantValue',true);

%% Compartment
comp_C = addcompartment(model,'V_c'      , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Central compartment (C) representing volume of the plasma ']);
comp_C = addcompartment(model,'V_p'      , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Peripheral compartment (P) representing volume of the peripheral tissues ']);
comp_C = addcompartment(model,'V_p_Oligo', 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Peripheral compartment (P) representing volume of the peripheral tissues ']);
comp_C = addcompartment(model,'V_sc'     , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Subcutaneous compartment (SC) representing volume of the SC space ']);
comp_C = addcompartment(model,'V_bsf'    , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Superficial Brain compartment (Bsf) representing the ISF volume of the brain accessible through diffusion from CSF ']);
comp_C = addcompartment(model,'V_bdp'    , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Deep Brain compartment (Bdp) representing the ISF volume of the brain not accessible through diffusion ']);
comp_C = addcompartment(model,'V_b'      , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Brain compartment (B) representing the brain homogenate with only derived numbers']);
comp_C = addcompartment(model,'V_csf'    , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['CSF compartment (CSF) representing volume of the CSF ']);
comp_C = addcompartment(model,'V_scord'  , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Spinal Cord compartment (SCord) representing the ISF volume of the spinal cord ']);
comp_C = addcompartment(model,'V_cv'     , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Cerebroventricular (CV) compartment representing volume of the CV space only used for ICV injections and assume anything not absorbed to leave into blood ']);
comp_C = addcompartment(model,'V_liver'  , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Liver (liver) compartment representing the isf volume of the liver']);
comp_C = addcompartment(model,'V_b_p'  , 1,'CapacityUnits','milliliter','ConstantCapacity',false); set(comp_C,'Notes',['Brain peripheral compartment']);
P = addparameter(model,'v_c'      , 1 ,'ValueUnits','milliliter','ConstantValue',true); 
P = addparameter(model,'v_p'      , 1 ,'ValueUnits','milliliter','ConstantValue',true); 
P = addparameter(model,'v_p_Oligo', 1 ,'ValueUnits','milliliter','ConstantValue',true); 
P = addparameter(model,'v_sc'     , 1 ,'ValueUnits','milliliter','ConstantValue',true); 
P = addparameter(model,'v_b'      , 1 ,'ValueUnits','milliliter','ConstantValue',true);  
P = addparameter(model,'v_csf'    , 1 ,'ValueUnits','milliliter','ConstantValue',true);  
P = addparameter(model,'v_scord'  , 1 ,'ValueUnits','milliliter','ConstantValue',true);  
P = addparameter(model,'v_cv'     , 1 ,'ValueUnits','milliliter','ConstantValue',true);  
P = addparameter(model,'v_liver'  , 1 ,'ValueUnits','milliliter','ConstantValue',true); 
P = addparameter(model,'v_b_p'  , 1 ,'ValueUnits','milliliter','ConstantValue',true);

P = addparameter(model,'f_sf2tot', 0.99 ,'ValueUnits','dimensionless','ConstantValue',true);  
rule = addrule(model,'V_c     = v_c'    ,'repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_p     = v_p'    ,'repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_sc    = v_sc'   ,'repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_csf   = v_csf'  ,'repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_scord = v_scord*f_isf2b','repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_cv    = v_cv'   ,'repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_b     = v_b'    ,'repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_p_Oligo = v_p_Oligo'       ,'repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_liver   = v_liver'         ,'repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_bsf     = v_b*   f_sf2tot*f_isf2b ','repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_bdp     = v_b*(1-f_sf2tot)*f_isf2b','repeatedAssignment'); set(rule,'Tag',' inVivo');
rule = addrule(model,'V_b_p = v_b_p','repeatedAssignment'); set(rule,'Tag',' inVivo');

display("Brain comp volumes represent isf. Check how your drug and target concentrations from brain homogenate fit the model baselines and outputs");

%% Species


% Oligo
S = addspecies(sbioselect(model,'name','V_c'),      'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_p_Oligo'),'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_sc'),     'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_liver'),  'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bsf'),    'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bdp'),    'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_b'),      'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_csf'),    'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_cv' ),    'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_scord' ), 'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_b_p' ), 'Oligo_nM'           ,0      ,'InitialAmountUnits','nanomolarity');
% Oligo in the cell
S = addspecies(sbioselect(model,'name','V_liver'),  'OligoCell_nM'       ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bsf'),    'OligoCell_nM'       ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bdp'),    'OligoCell_nM'       ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_b'),      'OligoCell_nM'       ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_scord'),  'OligoCell_nM'       ,0      ,'InitialAmountUnits','nanomolarity');
% mRNA
S = addspecies(sbioselect(model,'name','V_bsf'),    'mRNA_perc'      ,100      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bdp'),    'mRNA_perc'      ,100      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_liver'),  'mRNA_perc'      ,100      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_scord'),  'mRNA_perc'      ,100      ,'InitialAmountUnits','nanomolarity');
% protein
S = addspecies(sbioselect(model,'name','V_bsf'),    'Prot_perc'      ,100      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bdp'),    'Prot_perc'      ,100      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_liver'),  'Prot_perc'      ,100      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_scord'),  'Prot_perc'      ,100      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_csf'),    'Prot_perc'      ,100      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_csf'),    'Prot'      ,100      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_c'),      'Prot_perc'      ,100      ,'InitialAmountUnits','nanomolarity');

%% observables (derived from the primary species)

% Calculate AUCs
% S = addspecies(sbioselect(model,'name','V_c'),      'AUC_D_nMhr'       ,0      ,'InitialAmountUnits','nanomolarity');
% S = addspecies(sbioselect(model,'name','V_csf'),    'AUC_D_nMhr'       ,0      ,'InitialAmountUnits','nanomolarity');
% S = addspecies(sbioselect(model,'name','V_bdp'),    'AUC_D_nMhr'       ,0      ,'InitialAmountUnits','nanomolarity');
% S = addspecies(sbioselect(model,'name','V_bsf'),    'AUC_D_nMhr'       ,0      ,'InitialAmountUnits','nanomolarity');

% Drug derivatives
% Plasma

S = addspecies(sbioselect(model,'name','V_c'),    'Oligo_tot_nM'       ,0      ,'InitialAmountUnits','nanomolarity');

addrule(model,'V_c.Oligo_tot_nM  = V_c.Oligo_nM ','repeatedAssignment');
S = addspecies(sbioselect(model,'name','V_c'),    'Oligo_tot_ngml'         ,0      ,'InitialAmountUnits','nanogram/milliliter');
addrule(model,'V_c.Oligo_tot_ngml = V_c.Oligo_tot_nM*MW_D','repeatedAssignment');

% Superficial Brain
S = addspecies(sbioselect(model,'name','V_bsf'),    'Oligo_tot_nM'        ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bsf'),    'Oligo_tot_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bsf'),    'OligoCell_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
addrule(model,'V_bsf.Oligo_tot_nM  = V_bsf.Oligo_nM + V_bsf.OligoCell_nM','repeatedAssignment');
addrule(model,'V_bsf.Oligo_tot_homog_nM = V_bsf.Oligo_tot_nM*f_isf2b','repeatedAssignment');
addrule(model,'V_bsf.OligoCell_homog_nM = V_bsf.OligoCell_nM*f_isf2b','repeatedAssignment');
S = addspecies(sbioselect(model,'name','V_bsf'),    'Oligo_tot_homog_ngml '  ,0      ,'InitialAmountUnits','nanogram/milliliter');
addrule(model,'V_bsf.Oligo_tot_homog_ngml = V_bsf.Oligo_tot_homog_nM*MW_D'  ,'repeatedAssignment');
S = addspecies(sbioselect(model,'name','V_bsf'),    'Oligo_csf_norm '  ,0      ,'InitialAmountUnits','nanomolarity');
addrule(model,'V_bsf.Oligo_csf_norm = V_bsf.Oligo_tot_homog_nM/(V_csf.Oligo_tot_nM*one_nM)'  ,'repeatedAssignment');

% Deep Brain

S = addspecies(sbioselect(model,'name','V_bdp'),    'Oligo_tot_nM'        ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bdp'),    'Oligo_tot_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_bdp'),    'OligoCell_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
addrule(model,'V_bdp.Oligo_tot_nM  =  V_bdp.Oligo_nM + V_bdp.OligoCell_nM','repeatedAssignment');
addrule(model,'V_bdp.Oligo_tot_homog_nM = V_bdp.Oligo_tot_nM*f_isf2b','repeatedAssignment');
addrule(model,'V_bdp.OligoCell_homog_nM = V_bdp.OligoCell_nM*f_isf2b','repeatedAssignment');
S = addspecies(sbioselect(model,'name','V_bdp'),    'Oligo_tot_homog_ngml'  ,0      ,'InitialAmountUnits','nanogram/milliliter');
addrule(model,'V_bdp.Oligo_tot_homog_ngml = V_bdp.Oligo_tot_homog_nM*MW_D','repeatedAssignment');

% Whole Brain
addrule(model,'V_b.Oligo_nM     = (1-f_sf2tot)*V_bdp.Oligo_nM     + f_sf2tot*V_bsf.Oligo_nM'    ,'repeatedAssignment');
addrule(model,'V_b.OligoCell_nM = (1-f_sf2tot)*V_bdp.OligoCell_nM + f_sf2tot*V_bsf.OligoCell_nM','repeatedAssignment');
S = addspecies(sbioselect(model,'name','V_b'),    'Oligo_tot_nM'        ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_b'),    'Oligo_tot_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_b'),    'OligoCell_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
addrule(model,'V_b.Oligo_tot_nM       = (1-f_sf2tot)*V_bdp.Oligo_tot_nM       + f_sf2tot*V_bsf.Oligo_tot_nM'      ,'repeatedAssignment');
addrule(model,'V_b.Oligo_tot_homog_nM = (1-f_sf2tot)*V_bdp.Oligo_tot_homog_nM + f_sf2tot*V_bsf.Oligo_tot_homog_nM','repeatedAssignment');
addrule(model,'V_b.OligoCell_homog_nM = V_b.OligoCell_nM*f_isf2b','repeatedAssignment');
S = addspecies(sbioselect(model,'name','V_b'),    'Oligo_tot_homog_perc'  ,0      ,'InitialAmountUnits','nanomolarity');
addrule(model,'V_b.Oligo_tot_homog_perc = 100*V_b.Oligo_tot_homog_nM*v_b*MW_D/Dose_0'        ,'repeatedAssignment');
% Spinal cord
S = addspecies(sbioselect(model,'name','V_scord'),    'Oligo_tot_nM'        ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_scord'),    'Oligo_tot_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_scord'),    'OligoCell_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
addrule(model,'V_scord.Oligo_tot_nM  = V_scord.Oligo_nM + V_scord.OligoCell_nM','repeatedAssignment');
addrule(model,'V_scord.Oligo_tot_homog_nM = V_scord.Oligo_tot_nM*f_isf2b','repeatedAssignment');
addrule(model,'V_scord.OligoCell_homog_nM = V_scord.OligoCell_nM*f_isf2b','repeatedAssignment');
S = addspecies(sbioselect(model,'name','V_scord'),    'Oligo_tot_homog_ngml'  ,0      ,'InitialAmountUnits','nanogram/milliliter');
addrule(model,'V_scord.Oligo_tot_homog_ngml = V_scord.Oligo_tot_homog_nM*MW_D '  ,'repeatedAssignment');

% CSF
S = addspecies(sbioselect(model,'name','V_csf'),    'Oligo_tot_nM'       ,0      ,'InitialAmountUnits','nanomolarity');
addrule(model,'V_csf.Oligo_tot_nM  = V_csf.Oligo_nM ','repeatedAssignment');
S = addspecies(sbioselect(model,'name','V_csf'),    'Oligo_tot_ngml'         ,0      ,'InitialAmountUnits','nanogram/milliliter');
addrule(model,'V_csf.Oligo_tot_ngml = V_csf.Oligo_tot_nM*MW_D ','repeatedAssignment');
addrule(model,'V_csf.Prot_perc = 100*V_csf.Prot/(one_nM*Prot_perc_bl*kdiff_b2csf_Prot*(V_bsf+V_bdp)/(V_csf*(0.693/t_hl_csf)))','repeatedAssignment');
addrule(model,'V_csf.Prot = (Prot_perc_bl*kdiff_b2csf_Prot*(V_bsf+V_bdp)/(V_csf*(0.693/t_hl_csf)))','initialAssignment');

% Liver
S = addspecies(sbioselect(model,'name','V_liver'),    'Oligo_tot_nM'        ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_liver'),    'Oligo_tot_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
S = addspecies(sbioselect(model,'name','V_liver'),    'OligoCell_homog_nM'  ,0      ,'InitialAmountUnits','nanomolarity');
addrule(model,'V_liver.Oligo_tot_nM  = V_liver.Oligo_nM + V_liver.OligoCell_nM','repeatedAssignment');
addrule(model,'V_liver.Oligo_tot_homog_nM = V_liver.Oligo_tot_nM*f_isf2liver','repeatedAssignment');
addrule(model,'V_liver.OligoCell_homog_nM = V_liver.OligoCell_nM*f_isf2b','repeatedAssignment');


%% Parameters  

% Basics
P = addparameter(model,'NA',6.02214199e23,'ValueUnits','molecule/mole','ConstantValue',true); set(P,'Notes',['Avogadro''s number']);
P = addparameter(model,'BW_mouse'   , 0.025 ,'ValueUnits','kilogram','ConstantValue',false);  set(P,'Notes',['Body weight of mouse']); 
P = addparameter(model,'BW_monkey'  , 3.8   ,'ValueUnits','kilogram','ConstantValue',false);  set(P,'Notes',['Body weight of monkey']); 
P = addparameter(model,'BW_human'   ,70     ,'ValueUnits','kilogram','ConstantValue',false);  set(P,'Notes',['Body weight of human']); 

% flags
P = addparameter(model,'flag_mouse' ,   0 ,'ValueUnits','dimensionless','ConstantValue',true);
P = addparameter(model,'flag_monkey',   0 ,'ValueUnits','dimensionless','ConstantValue',true);
P = addparameter(model,'flag_human' ,   0 ,'ValueUnits','dimensionless','ConstantValue',true);
% flags for molecules
P = addparameter(model,'flag_HLE'    ,   0 ,'ValueUnits','dimensionless','ConstantValue',true);

% MW of the drug
P = addparameter(model,'MW_Oligo'  , 14 ,'ValueUnits','kilogram/mole','ConstantValue',false);  set(P,'Notes',['molecular weight of the Oligo']);
P = addparameter(model,'MW_D'      ,110 ,'ValueUnits','kilogram/mole','ConstantValue',false);  set(P,'Notes',['molecular weight of the molecule']); 
addrule(model,'MW_D = MW_Oligo ','repeatedAssignment');

% PK Parameters
P = addparameter(model,'CL_mouse'      ,   10 ,'ValueUnits','milliliter/day','ConstantValue',true);
P = addparameter(model,'CL_monkey'     ,   10 ,'ValueUnits','milliliter/day','ConstantValue',true);
P = addparameter(model,'CL_human'      ,   10 ,'ValueUnits','milliliter/day','ConstantValue',true);
P = addparameter(model,'eta_CL'        ,    1 ,'ValueUnits','dimensionless' ,'ConstantValue',true);
P = addparameter(model,'nAllo_CL_m2h'  , 0.90 ,'ValueUnits','dimensionless' ,'ConstantValue',true);
P = addparameter(model,'nAllo_CL_c2h'  , 0.85 ,'ValueUnits','dimensionless' ,'ConstantValue',true);
P = addparameter(model,'flag_CL_m2h'   , 0    ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'flag_CL_c2h'   , 0    ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'flag_CL_h2h'   , 0    ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'CLd'     ,   10 ,'ValueUnits','milliliter/day' ,'ConstantValue',true);
P = addparameter(model,'f_HLE'   ,    3 ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'f_Rec'   ,    1 ,'ValueUnits','dimensionless'  ,'ConstantValue',true); % not in the equations yet
P = addparameter(model,'Vm_TfR'  ,    2 ,'ValueUnits','nanomole/hour'  ,'ConstantValue',true);
P = addparameter(model,'Km_TfR'  ,   10 ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);
P = addparameter(model,'f_rec'   ,   1  ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'ka'      , 0.01 ,'ValueUnits','1/hour'         ,'ConstantValue',true);
P = addparameter(model,'fa'      ,   1  ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
% oligo PK, stability, and conjugation related parameters
P = addparameter(model,'CL_Oligo'            ,   1000      ,'ValueUnits','milliliter/day' ,'ConstantValue',false);
P = addparameter(model,'f_CL_Oligo_2liver'   ,    0.2      ,'ValueUnits','dimensionless'  ,'ConstantValue',false);
P = addparameter(model,'CLd_Oligo'           ,   1000      ,'ValueUnits','milliliter/day' ,'ConstantValue',false);
P = addparameter(model,'t_hl_deconj'         ,    168      ,'ValueUnits','hour'           ,'ConstantValue',true);
P = addparameter(model,'t_hl_breakdown'      ,    72       ,'ValueUnits','hour'           ,'ConstantValue',true);
P = addparameter(model,'t_hl_Oligo_liver'    ,    7        ,'ValueUnits','day'            ,'ConstantValue',true);
P = addparameter(model,'t_hl_Oligo_b'        ,    90       ,'ValueUnits','day'            ,'ConstantValue',true);


% Brain/CSF uptake Parameters
P = addparameter(model,'Kp_b2c_D'         ,    0.0001  ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'Kp_csf2c_D'       ,    0.001   ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'kdeg_D_b'         ,    0.0270  ,'ValueUnits','1/hour'         ,'ConstantValue',true);
P = addparameter(model,'f_TfR_BBB'        ,    0.03    ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'f_TfR_scord'      ,    0.001   ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'kdeg_TfR_D_b'     ,    0.1540  ,'ValueUnits','1/hour'         ,'ConstantValue',true);
P = addparameter(model,'Km_kdeg_D_b'      ,    5500    ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);
P = addparameter(model,'f_isf2b'          ,    0.2     ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'f_isf2liver'      ,    0.2     ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'t_hl_csf'         ,    8       ,'ValueUnits','hour'           ,'ConstantValue',true);
P = addparameter(model,'kdiff_csf2bsf_D'  ,    0.075*1 ,'ValueUnits','milliliter/day'          ,'ConstantValue',true);
P = addparameter(model,'kdiff_csf2bdp_D'  ,    0.012*1,'ValueUnits','milliliter/day'          ,'ConstantValue',true);
P = addparameter(model,'kdiff_csf2b_p_D'  ,    0.075*1 ,'ValueUnits','milliliter/day'          ,'ConstantValue',true);
P = addparameter(model,'kdiff_csf2scord_D',    0.075/10 ,'ValueUnits','milliliter/day'          ,'ConstantValue',true);
P = addparameter(model,'ka_ICV'           ,    1       ,'ValueUnits','1/day'          ,'ConstantValue',true);
P = addparameter(model,'fa_ICV_bsf'       ,    0  ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'fa_ICV_bdp'       ,    0  ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'fa_ICV'           ,    0       ,'ValueUnits','dimensionless'  ,'ConstantValue',false);
P = addparameter(model,'f_cellUptake'     ,    0.35    ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'t_hl_cellUptake'  ,    2       ,'ValueUnits','day'            ,'ConstantValue',true);
addrule(model,'fa_ICV = fa_ICV_bsf*f_sf2tot + fa_ICV_bdp*(1-f_sf2tot)','repeatedAssignment');

% mRNA 
P = addparameter(model,'t_hl_mRNA'          ,    20      ,'ValueUnits','hour'           ,'ConstantValue',true);set(P,'Notes',['20 hr for tau (https://pubmed.ncbi.nlm.nih.gov/10527457/)']);
P = addparameter(model,'IC50_mRNA_Oligo_b'  ,    50      ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);
P = addparameter(model,'IC50_mRNA_Oligo_l'  ,    50      ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);
P = addparameter(model,'n_mRNA'             ,    1       ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'mRNA_perc_bl'            ,    100     ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);

% Protein
P = addparameter(model,'t_hl_Prot_b'       ,    23      ,'ValueUnits','day'            ,'ConstantValue',true);set(P,'Notes',['20 hr for tau (doi:10.1016/j.neuron.2018.02.015)']);
P = addparameter(model,'t_hl_delay'       ,    23      ,'ValueUnits','day'            ,'ConstantValue',true);
P = addparameter(model,'t_hl_Prot_scord'   ,    23      ,'ValueUnits','day'            ,'ConstantValue',true);set(P,'Notes',['20 hr for tau (doi:10.1016/j.neuron.2018.02.015)']);
P = addparameter(model,'t_hl_Prot_l'       ,    23      ,'ValueUnits','day'            ,'ConstantValue',true);set(P,'Notes',['20 hr for tau (doi:10.1016/j.neuron.2018.02.015)']);
P = addparameter(model,'t_hl_Prot_c'       ,    0.4     ,'ValueUnits','day'            ,'ConstantValue',true);set(P,'Notes',['20 hr for tau (doi:10.1016/j.neuron.2018.02.015)']);
P = addparameter(model,'n_Prot'            ,    1       ,'ValueUnits','dimensionless'  ,'ConstantValue',true);
P = addparameter(model,'Prot_perc_bl'      ,    100     ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);
P = addparameter(model,'Prot_bl_c'         ,    10      ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);
P = addparameter(model,'Prot_bl_csf'       ,    100     ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);
P = addparameter(model,'Prot_bl_bsf'       ,    100000  ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);
P = addparameter(model,'Prot_bl_bdp'       ,    100000  ,'ValueUnits','nanomolarity'   ,'ConstantValue',true);
P = addparameter(model,'kdiff_b2csf_Prot'       ,    3e-6  ,'ValueUnits','1/day'   ,'ConstantValue',true);


% Dummy
P = addparameter(model,'one_nM'       ,    1  ,'ValueUnits','1/nanomolarity'   ,'ConstantValue',true);

% 
% Tau levels and tau content varies by compartment and disease state. Average levels of tau in the brain are 175,000–220,000 times higher than levels of tau in the CSF. Total CSF tau is about 10 times higher than total plasma tau.
% add 

% Dose 
P = addparameter(model,'Dose_0'       ,    1  ,'ValueUnits','milligram/nanomolarity'   ,'ConstantValue',true); % for converting conc to % injected dose

%% Reactions
% SC absorption
reaction = addreaction(model,'V_sc.Oligo_nM    -> V_c.Oligo_nM');                 set(reaction,'ReactionRate','ka*   fa *V_sc.Oligo_nM');
reaction = addreaction(model,'V_sc.Oligo_nM    -> null');                         set(reaction,'ReactionRate','ka*(1-fa)*V_sc.Oligo_nM');

% Systemic PK
reaction = addreaction(model,'V_c.Oligo_nM -> null');                 set(reaction,'ReactionRate','CL_Oligo*V_c.Oligo_nM/V_c');
reaction = addreaction(model,'V_c.Oligo_nM -> V_p_Oligo.Oligo_nM');   set(reaction,'ReactionRate','CLd_Oligo*(V_c.Oligo_nM - V_p_Oligo.Oligo_nM)');


% Liver PK
reaction = addreaction(model,'null -> V_liver.Oligo_nM');              set(reaction,'ReactionRate','f_CL_Oligo_2liver*CL_Oligo*V_c.Oligo_nM/V_c');
% From Oligo to cellular fraction that is active (capture delay in mRNA onset)
reaction = addreaction(model,'V_liver.Oligo_nM     -> V_liver.OligoCell_nM');  set(reaction,'ReactionRate','(0.693/t_hl_cellUptake)*V_liver.Oligo_nM');
reaction = addreaction(model,'V_liver.OligoCell_nM -> null');                  set(reaction,'ReactionRate','(0.693/t_hl_Oligo_liver)*V_liver.OligoCell_nM');
% TASK: Add TfR effect in distriution in liver hepatocyte vs other cells


% Brain Uptake PK
    % UPTAKE
    % non-specific diffusion with different rates to bsf vs bdp (two-way)
    reaction = addreaction(model,'V_csf.Oligo_nM      -> V_bsf.Oligo_nM');      set(reaction,'ReactionRate','kdiff_csf2bsf_D*(V_csf.Oligo_nM-V_bsf.Oligo_nM)' );
    reaction = addreaction(model,'V_csf.Oligo_nM      -> V_bdp.Oligo_nM');      set(reaction,'ReactionRate','kdiff_csf2bdp_D*(V_csf.Oligo_nM-V_bdp.Oligo_nM)');
    reaction = addreaction(model,'V_csf.Oligo_nM      -> V_b_p.Oligo_nM');      set(reaction,'ReactionRate','kdiff_csf2b_p_D*(V_csf.Oligo_nM-V_b_p.Oligo_nM)' );
    
    
        
    % DEGRADATION
    % non-specific degradation in both bsf and bdp (TASK re-check the oligo ns deg)
    reaction = addreaction(model,'V_bsf.Oligo_nM -> null');       set(reaction,'ReactionRate','(0.693/t_hl_Oligo_b)*V_bsf.Oligo_nM');
    reaction = addreaction(model,'V_bdp.Oligo_nM -> null');       set(reaction,'ReactionRate','(0.693/t_hl_Oligo_b)*V_bdp.Oligo_nM');
    
    % CELLULAR UPTAKE
    % From Oligo to cellular fraction that is active (capture delay in mRNA onset and differentiate cellular hl)
    reaction = addreaction(model,'V_bsf.Oligo_nM     -> V_bsf.OligoCell_nM');    set(reaction,'ReactionRate','(0.693/t_hl_cellUptake)*V_bsf.Oligo_nM');
    reaction = addreaction(model,'V_bdp.Oligo_nM     -> V_bdp.OligoCell_nM');    set(reaction,'ReactionRate','(0.693/t_hl_cellUptake)*V_bdp.Oligo_nM');
    reaction = addreaction(model,'V_bsf.OligoCell_nM -> null');                  set(reaction,'ReactionRate','(0.693/t_hl_Oligo_b)*V_bsf.OligoCell_nM');
    reaction = addreaction(model,'V_bdp.OligoCell_nM -> null');                  set(reaction,'ReactionRate','(0.693/t_hl_Oligo_b)*V_bdp.OligoCell_nM');

    % ICV absorption
    reaction = addreaction(model,'V_cv.Oligo_nM    -> V_bsf.Oligo_nM');     set(reaction,'ReactionRate','ka_ICV*   fa_ICV_bsf*   f_sf2tot                      *V_cv.Oligo_nM');
    reaction = addreaction(model,'V_cv.Oligo_nM    -> V_bdp.Oligo_nM');     set(reaction,'ReactionRate','ka_ICV*   fa_ICV_bdp*(1-f_sf2tot)                     *V_cv.Oligo_nM');
    reaction = addreaction(model,'V_cv.Oligo_nM    -> V_c.Oligo_nM');       set(reaction,'ReactionRate','ka_ICV*(1-fa_ICV_bsf*f_sf2tot-fa_ICV_bdp*(1-f_sf2tot))*V_cv.Oligo_nM');
    
% Spinal Cord V_scord
  % UPTAKE
    % non-specific diffusion with different rates to bsf vs bdp
    reaction = addreaction(model,'V_csf.Oligo_nM      -> V_scord.Oligo_nM');      set(reaction,'ReactionRate','kdiff_csf2scord_D*(V_csf.Oligo_nM-V_scord.Oligo_nM)');
    % DEGRADATION
    % non-specific degradation in spinal cord
    reaction = addreaction(model,'V_scord.Oligo_nM -> null');       set(reaction,'ReactionRate','(0.693/t_hl_Oligo_b)*V_scord.Oligo_nM');
    % CELLULAR UPTAKE
    % From Oligo to cellular fraction that is active (capture delay in mRNA onset)
    reaction = addreaction(model,'V_scord.Oligo_nM     -> V_scord.OligoCell_nM');    set(reaction,'ReactionRate','(0.693/t_hl_cellUptake)*V_scord.Oligo_nM');
    reaction = addreaction(model,'V_scord.OligoCell_nM -> null');                    set(reaction,'ReactionRate','(0.693/t_hl_Oligo_b)*V_scord.OligoCell_nM');


% CSF Uptake PK
reaction = addreaction(model,'V_c.Oligo_nM      -> V_csf.Oligo_nM');      set(reaction,'ReactionRate','Kp_csf2c_D*0.693/t_hl_csf*V_csf/V_c*V_c.Oligo_nM');
reaction = addreaction(model,'V_csf.Oligo_nM    -> V_c.Oligo_nM ');       set(reaction,'ReactionRate','           0.693/t_hl_csf*V_csf.Oligo_nM');

% Brain PD
reaction = addreaction(model,'null -> V_bsf.mRNA_perc');    set(reaction,'ReactionRate','(0.693/t_hl_mRNA)*mRNA_perc_bl/(1+(V_bsf.OligoCell_homog_nM/IC50_mRNA_Oligo_b)^n_mRNA)   -(0.693/t_hl_mRNA)*V_bsf.mRNA_perc');
reaction = addreaction(model,'null -> V_bsf.Prot_perc');    set(reaction,'ReactionRate','(0.693/t_hl_Prot_b)*Prot_perc_bl*(V_bsf.mRNA_perc /mRNA_perc_bl)^n_Prot                   -(0.693/t_hl_Prot_b)*V_bsf.Prot_perc');
reaction = addreaction(model,'null -> V_bdp.mRNA_perc');    set(reaction,'ReactionRate','(0.693/t_hl_mRNA)*mRNA_perc_bl/(1+(V_bdp.OligoCell_homog_nM/IC50_mRNA_Oligo_b)^n_mRNA)   -(0.693/t_hl_mRNA)*V_bdp.mRNA_perc');
reaction = addreaction(model,'null -> V_bdp.Prot_perc');    set(reaction,'ReactionRate','(0.693/t_hl_Prot_b)*Prot_perc_bl*(V_bdp.mRNA_perc /mRNA_perc_bl)^n_Prot                   -(0.693/t_hl_Prot_b)*V_bdp.Prot_perc');
% Liver PD
reaction = addreaction(model,'null -> V_liver.mRNA_perc');    set(reaction,'ReactionRate','(0.693/t_hl_mRNA)*mRNA_perc_bl/(1+(V_liver.OligoCell_homog_nM/IC50_mRNA_Oligo_l)^n_mRNA)   -(0.693/t_hl_mRNA)*V_liver.mRNA_perc');
reaction = addreaction(model,'null -> V_liver.Prot_perc');    set(reaction,'ReactionRate','(0.693/t_hl_Prot_l)*Prot_perc_bl*(V_liver.mRNA_perc /mRNA_perc_bl)^n_Prot                   -(0.693/t_hl_Prot_l)*V_liver.Prot_perc');
% Spinal Cord PD
reaction = addreaction(model,'null -> V_scord.mRNA_perc');    set(reaction,'ReactionRate','(0.693/t_hl_mRNA)*mRNA_perc_bl/(1+(V_scord.OligoCell_homog_nM/IC50_mRNA_Oligo_b)^n_mRNA)   -(0.693/t_hl_mRNA)*V_scord.mRNA_perc');
reaction = addreaction(model,'null -> V_scord.Prot_perc');    set(reaction,'ReactionRate','(0.693/t_hl_Prot_scord)*Prot_perc_bl*(V_scord.mRNA_perc /mRNA_perc_bl)^n_Prot                   -(0.693/t_hl_Prot_scord)*V_scord.Prot_perc');

% CSF Protein dynamics
reaction = addreaction(model,'V_bsf.Prot_perc -> V_csf.Prot');    set(reaction,'ReactionRate','kdiff_b2csf_Prot*V_bsf.Prot_perc');
reaction = addreaction(model,'V_bdp.Prot_perc -> V_csf.Prot');    set(reaction,'ReactionRate','kdiff_b2csf_Prot*V_bdp.Prot_perc');
reaction = addreaction(model,'V_csf.Prot -> null');    set(reaction,'ReactionRate','0.693/t_hl_csf*V_csf.Prot');



% Protein dynamics
% 1. copy protein_perc but with actual baseline
% 2. add brain to CSF (use similar to PK diffusion out + volumes in contribution to CSF PD)
% 3. add CSF to plasma uptake flux based on the size of the protein
% 4. add CSF turnover clearing of protein from CSF
% 5. add clearance from plasma to maintain plasma baseline

%% Dose Parameters  
% Drug dose IV
P = addparameter(model,'dose_mgkg_IV'           ,0  ,'ValueUnits','milligram'    ,'ConstantValue',true); 
P = addparameter(model,'dose_Oligo_nmole_IV'    ,0  ,'ValueUnits','picomole'     ,'ConstantValue',false);  
addrule(model,'dose_Oligo_nmole_IV    = dose_mgkg_IV / MW_D','initialAssignment');
P = addparameter(model,'startTime_day_IV'  ,0  ,'ValueUnits','day'          ,'ConstantValue',true); 
P = addparameter(model,'interval_day_IV'   ,0  ,'ValueUnits','day'          ,'ConstantValue',true); 
P = addparameter(model,'reapeatCount_IV'   ,0  ,'ValueUnits','dimensionless','ConstantValue',true); 
P = addparameter(model,'infusionDur_min_IV',0.1,'ValueUnits','minute'       ,'ConstantValue',true); 

% Drug dose SC
P = addparameter(model,'dose_mgkg_SC'           ,0  ,'ValueUnits','milligram'    ,'ConstantValue',true); 
P = addparameter(model,'dose_Oligo_nmole_SC'    ,0  ,'ValueUnits','picomole'     ,'ConstantValue',false);
addrule(model,'dose_Oligo_nmole_SC    = dose_mgkg_SC / MW_D','initialAssignment');
P = addparameter(model,'startTime_day_SC'  ,0  ,'ValueUnits','day'          ,'ConstantValue',true); 
P = addparameter(model,'interval_day_SC'   ,0  ,'ValueUnits','day'          ,'ConstantValue',true); 
P = addparameter(model,'reapeatCount_SC'   ,0  ,'ValueUnits','dimensionless','ConstantValue',true); 
P = addparameter(model,'infusionDur_min_SC',0.1,'ValueUnits','minute'       ,'ConstantValue',true); 

% Drug dose IT
P = addparameter(model,'dose_mgkg_IT'           ,0  ,'ValueUnits','milligram'    ,'ConstantValue',true); 
P = addparameter(model,'dose_Oligo_nmole_IT'    ,0  ,'ValueUnits','picomole'     ,'ConstantValue',false); 
addrule(model,'dose_Oligo_nmole_IT    = dose_mgkg_IT / MW_D','initialAssignment');
P = addparameter(model,'startTime_day_IT'  ,0  ,'ValueUnits','day'          ,'ConstantValue',true); 
P = addparameter(model,'interval_day_IT'   ,0  ,'ValueUnits','day'          ,'ConstantValue',true); 
P = addparameter(model,'reapeatCount_IT'   ,0  ,'ValueUnits','dimensionless','ConstantValue',true); 
P = addparameter(model,'infusionDur_min_IT',0.1,'ValueUnits','minute'       ,'ConstantValue',true); 

% Drug dose ICV
P = addparameter(model,'dose_mgkg_ICV'           ,0  ,'ValueUnits','milligram'    ,'ConstantValue',true); 
P = addparameter(model,'dose_Oligo_nmole_ICV'    ,0  ,'ValueUnits','picomole'     ,'ConstantValue',false);
addrule(model,'dose_Oligo_nmole_ICV    = dose_mgkg_ICV / MW_D','initialAssignment');
P = addparameter(model,'startTime_day_ICV'  ,0  ,'ValueUnits','day'          ,'ConstantValue',true); 
P = addparameter(model,'interval_day_ICV'   ,0  ,'ValueUnits','day'          ,'ConstantValue',true); 
P = addparameter(model,'reapeatCount_ICV'   ,0  ,'ValueUnits','dimensionless','ConstantValue',true); 
P = addparameter(model,'infusionDur_min_ICV',0.1,'ValueUnits','minute'       ,'ConstantValue',true); 


%% diagnostics
display("Diagnostics:");
unitlessObjects = sbioselect(model, 'Units', '')


