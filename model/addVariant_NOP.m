function model = addVariant_NOP(model,data,args,options)
% addVariant_NOP  Define the SimBiology variants used by the NOP scenarios.
%   model = addVariant_NOP(model,data,args) adds the parameter variants applied
%   when simulating different molecules, species and study contexts: species
%   flag variants, per-molecule physiology variants (mouse/cyno/human), PK
%   variants (diffusion + peripheral volume), PD variants and IC50 variants.
%
%   Inputs:
%     model   - SimBiology model to add variants to.
%     data    - dataset table (study/group metadata).
%     args    - struct of project settings.
%     options - name-value options (Name).
%   Outputs:
%     model   - model with the variants added.

arguments
    model SimBiology.Model
    data table
    args struct
    options.Name string = ''
end

%% create a null variant ro deal with table storage of variants
variantObj = addvariant(model,'varNull');

% TASK: test if fitting the CL
%% create the species variant
variantObj = addvariant(model,'varSpecies_mouse');
addcontent(variantObj, {'parameter', 'flag_mouse'   , 'Value', 1});
addcontent(variantObj, {'parameter', 'flag_monkey'  , 'Value', 0});
addcontent(variantObj, {'parameter', 'flag_human'   , 'Value', 0});

% Mouse Brain volume - 0.485 ml
% Mouse csf volume - 35 ul % Spine - 14 ul (csf to spine ratio from cyno)
% Plasma - 0.944 ml
variantObj = addvariant(model,'varMol_Roche_mouse'); % 25 gm
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 6});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 19.4});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.4});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 1.8});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 406.317});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4});
addcontent(variantObj, {'parameter', 'ka_ICV', 'Value', 9.24}); % Calculated based on csf turnover
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.99});
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 30.13});
% ICV
addcontent(variantObj, {'parameter', 'fa_ICV_bdp', 'Value', 0});
addcontent(variantObj, {'parameter', 'fa_ICV_bsf', 'Value', 0.025});


variantObj = addvariant(model,'varMol_Biogen_mouse'); % 25 gm
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 19.4});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.4});
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 0.56});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 1.8});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4});
addcontent(variantObj, {'parameter', 'ka_ICV', 'Value', 9.24}); % Calculated based on csf turnover
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.99});
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 33.1309});
% ICV
addcontent(variantObj, {'parameter', 'fa_ICV_bdp', 'Value', 0});
addcontent(variantObj, {'parameter', 'fa_ICV_bsf', 'Value', 0.0201});


variantObj = addvariant(model,'varMol_RO7234292_mouse'); % 25 gm
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 19.4});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.4});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 1.8});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 182.8571});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4});
addcontent(variantObj, {'parameter', 'ka_ICV', 'Value', 9.24}); % Calculated based on csf turnover
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.99});
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 29.124}); % 28.8226
% ICV
addcontent(variantObj, {'parameter', 'fa_ICV_bdp', 'Value', 0});
addcontent(variantObj, {'parameter', 'fa_ICV_bsf', 'Value', 0.0737});

variantObj = addvariant(model,'varMol_Jacifusen_mouse'); % 25 gm
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 19.4});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.4});
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 0.56});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 1.8});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 300.67});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4});
addcontent(variantObj, {'parameter', 'ka_ICV', 'Value', 9.24}); % Calculated based on csf turnover
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.99});
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 57.6025});
% ICV
addcontent(variantObj, {'parameter', 'fa_ICV_bdp', 'Value', 0});
addcontent(variantObj, {'parameter', 'fa_ICV_bsf', 'Value', 0.0091});


% Cyno Brain volume - 94 ml, Cortex - 54 ml
% Cyno CSF volume - 15 ml
% Cyno spine volume - 6 ml
variantObj = addvariant(model,'varMol_BIIB080_cyno_Easton'); % 6.75 kg
addcontent(variantObj, {'parameter', 'MW_D'        , 'Value', 6});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 56});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 13.92});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.92});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 5.2});
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 1.5});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-3});
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.57});

variantObj = addvariant(model,'varMol_BIIB080_cyno_Monine'); % 3.8 kg
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 58.95});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 24.74});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 3.95});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 5.2});
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 1.58});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4}); % 1e-3
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.57}); % 0.7625
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 57.28}); %3214.5
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 36.22});

variantObj = addvariant(model,'varMol_cyno_BrownKM2021'); % 3.8 kg
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 13.5});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 58.95});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 24.74});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 3.95});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 5.2});
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 1.58});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4}); % 1e-3
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.57}); % 0.7625
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 161});
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 118});

variantObj = addvariant(model,'varMol_Tominersen_cyno'); % 3.8 kg
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 58.95});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 24.74});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 3.95});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 5.2});
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 1.58});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4}); % 1e-3
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.57}); % 0.7625
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 675.0230}); % 324
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 18.52});

variantObj = addvariant(model,'varMol_Jacifusen_cyno'); % 3.8 kg
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 58.95});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 24.74});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 3.95});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 5.2});
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 1.58});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4}); % 1e-3
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.57}); % 0.7625
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 2051.1}); % 324
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 57.6}); % 76.19

% Human Brain volume - 1450 ml
% Human csf volume - 100 ml 
% Human spine volume - 40 ml (ratio from cyno) 
% Human plasma volume - 3126 ml (Shah)
variantObj = addvariant(model,'varMol_BIIB080_human');
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 44.66});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 20.71});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.43});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 4.8}); % 6.98
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 0.57});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4});
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.8}); % 0.57
addcontent(variantObj, {'parameter', 'Dose_0', 'Value', 0.14}); % 0.14 for PET, % 1.7143
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 68.42});
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 36.22});

variantObj = addvariant(model,'varMol_human_BrownKM2021');
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 13.5});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 44.66});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 20.71});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.43});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 4.8}); % 6.98
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 0.57});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4});
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.8}); % 0.57
addcontent(variantObj, {'parameter', 'Dose_0', 'Value', 0.14}); % 0.14 for PET, % 1.7143
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 148.2});
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 118});

variantObj = addvariant(model,'varMol_Tominersen_human');
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 44.66});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 20.71});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.43});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 4.8}); % 6.98
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 0.57});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4});
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.8}); % 0.57
addcontent(variantObj, {'parameter', 'Dose_0', 'Value', 0.14}); % 0.14 for PET, % 1.7143
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 668.39});
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 24.161});


variantObj = addvariant(model,'varMol_Jacifusen_human');
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 44.66});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 20.71});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.43});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 4.8}); % 6.98
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 0.57});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4});
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.8}); % 0.57
addcontent(variantObj, {'parameter', 'Dose_0', 'Value', 0.14}); % 0.14 for PET, % 1.7143
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 2501.1});
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 57.6});

variantObj = addvariant(model,'varMol_get_dose');
addcontent(variantObj, {'parameter', 'MW_Oligo'        , 'Value', 7});
addcontent(variantObj, {'parameter', 'v_c', 'Value', 44.66});
addcontent(variantObj, {'parameter', 'v_b', 'Value', 20.71});
addcontent(variantObj, {'parameter', 'v_csf', 'Value', 1.43});
addcontent(variantObj, {'parameter', 't_hl_csf', 'Value', 4.8}); % 6.98
addcontent(variantObj, {'parameter', 'v_scord', 'Value', 0.57});
addcontent(variantObj, {'parameter', 'Kp_csf2c_D', 'Value', 1e-4});
addcontent(variantObj, {'parameter', 'f_sf2tot', 'Value', 0.8}); % 0.57
addcontent(variantObj, {'parameter', 'Dose_0', 'Value', 0.14}); % 0.14 for PET, % 1.7143
addcontent(variantObj, {'parameter', 't_hl_Oligo_b', 'Value', 36.22});


%% PK Variants
%% Mouse (kdiff + v_b_p)
variantObj = addvariant(model,'mouse_kdiff_and_v_b_p');
addcontent(variantObj, {'parameter', 'kdiff_csf2bsf_D', 'Value', 0.99});
addcontent(variantObj, {'parameter', 'kdiff_csf2bdp_D', 'Value', 0.006});
addcontent(variantObj, {'parameter', 'kdiff_csf2scord_D', 'Value', 0.05});
addcontent(variantObj, {'parameter', 'kdiff_csf2b_p_D', 'Value', 0.04});
addcontent(variantObj, {'parameter', 'v_b_p', 'Value', 1.99});
% 2 compartment parameters
addcontent(variantObj, {'parameter', 'v_p_Oligo', 'Value', 1000});
addcontent(variantObj, {'parameter', 'CLd_Oligo', 'Value', 27});
addcontent(variantObj, {'parameter', 'CL_Oligo', 'Value', 1000});
% Brain
addcontent(variantObj, {'parameter', 'f_isf2b', 'Value', 0.2});
addcontent(variantObj, {'parameter', 't_hl_cellUptake', 'Value', 1});

%% Cyno (kdiff + v_b_p)
variantObj = addvariant(model,'cyno_kdiff_and_v_b_p');
addcontent(variantObj, {'parameter', 'kdiff_csf2bsf_D', 'Value', 0.73});
addcontent(variantObj, {'parameter', 'kdiff_csf2bdp_D', 'Value', 0.32});
addcontent(variantObj, {'parameter', 'kdiff_csf2scord_D', 'Value', 0.14});
addcontent(variantObj, {'parameter', 'kdiff_csf2b_p_D', 'Value', 0.04});
addcontent(variantObj, {'parameter', 'v_b_p', 'Value', 1.99});
% 2 compartment parameters
addcontent(variantObj, {'parameter', 'v_p_Oligo', 'Value', 1000});
addcontent(variantObj, {'parameter', 'CLd_Oligo', 'Value', 43.032});
addcontent(variantObj, {'parameter', 'CL_Oligo', 'Value', 1000});
% Brain
addcontent(variantObj, {'parameter', 'f_isf2b', 'Value', 0.2});
addcontent(variantObj, {'parameter', 't_hl_cellUptake', 'Value', 1});

%% Human (kdiff + v_b_p)
variantObj = addvariant(model,'human_kdiff_and_v_b_p');
addcontent(variantObj, {'parameter', 'kdiff_csf2bsf_D', 'Value', 0.63});
addcontent(variantObj, {'parameter', 'kdiff_csf2bdp_D', 'Value', 0.18});
addcontent(variantObj, {'parameter', 'kdiff_csf2scord_D', 'Value', 0.06});
addcontent(variantObj, {'parameter', 'kdiff_csf2b_p_D', 'Value', 0.008});
addcontent(variantObj, {'parameter', 'v_b_p', 'Value', 1});
% 2 compartment parameters
addcontent(variantObj, {'parameter', 'v_p_Oligo', 'Value', 1e4});
addcontent(variantObj, {'parameter', 'CLd_Oligo', 'Value', 2.09e3});
addcontent(variantObj, {'parameter', 'CL_Oligo', 'Value', 200});
% Brain
addcontent(variantObj, {'parameter', 'f_isf2b', 'Value', 0.2});
addcontent(variantObj, {'parameter', 't_hl_cellUptake', 'Value', 1});

%% Tominersen_Cyno (kdiff + v_b_p)
variantObj = addvariant(model,'Tominersen_cyno_kdiff_and_v_b_p');
addcontent(variantObj, {'parameter', 'kdiff_csf2bsf_D', 'Value', 0.73});
addcontent(variantObj, {'parameter', 'kdiff_csf2bdp_D', 'Value', 0.22});
addcontent(variantObj, {'parameter', 'kdiff_csf2scord_D', 'Value', 0.14});
addcontent(variantObj, {'parameter', 'kdiff_csf2b_p_D', 'Value', 0.04});
addcontent(variantObj, {'parameter', 'v_b_p', 'Value', 2.19});
% 2 compartment parameters
addcontent(variantObj, {'parameter', 'v_p_Oligo', 'Value', 1000});
addcontent(variantObj, {'parameter', 'CLd_Oligo', 'Value', 43.032});
addcontent(variantObj, {'parameter', 'CL_Oligo', 'Value', 1000});
% Brain
addcontent(variantObj, {'parameter', 'f_isf2b', 'Value', 0.2});
addcontent(variantObj, {'parameter', 't_hl_cellUptake', 'Value', 1});

%% PD Variants
%% BIIB080 Mouse PD
variantObj = addvariant(model,'BIIB080_Mouse_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 10});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 117.23});

%% BIIB080 Cyno PD
variantObj = addvariant(model,'BIIB080_Cyno_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 16});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 1628.4});

%% BIIB080 Human PD
variantObj = addvariant(model,'BIIB080_Human_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 114.44});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 162.68});

%% APP siRNA Cyno PD
variantObj = addvariant(model,'APP_siRNA_Cyno_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 4.15});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 0.58});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 161});

%% APP siRNA Human PD
variantObj = addvariant(model,'APP_siRNA_Human_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 3.09});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 113.75});

%% Tominersen Mouse PD
variantObj = addvariant(model,'Tominersen_Mouse_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 0.3145});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 1698.5});

%% Tominersen Cyno PD
variantObj = addvariant(model,'Tominersen_Cyno_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 0.75});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1.41});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 675.02});

%% Tominersen Human PD
variantObj = addvariant(model,'Tominersen_Human_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 50}); %50
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1}); %1
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 826.34}); % 826.34

%% Jacifusen Mouse PD
variantObj = addvariant(model,'Jacifusen_Mouse_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 300.67});

%% Jacifusen Cyno PD
variantObj = addvariant(model,'Jacifusen_Cyno_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 1683});

%% Jacifusen Human PD
variantObj = addvariant(model,'Jacifusen_Human_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 52.5});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 281.46});

%% Human PD get dose
variantObj = addvariant(model,'get_dose_Human_PD');
addcontent(variantObj, {'parameter', 't_hl_mRNA', 'Value', 24});
addcontent(variantObj, {'parameter', 't_hl_Prot_b', 'Value', 114.44});
addcontent(variantObj, {'parameter', 'n_mRNA', 'Value', 1});
addcontent(variantObj, {'parameter', 'n_Prot', 'Value', 1});
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 117.23});

%% IC50 Variants
%% BIIB080 Mouse IC50
variantObj = addvariant(model,'BIIB080_Mouse_IC50');
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 117.23});

%% BIIB080 Cyno IC50
variantObj = addvariant(model,'BIIB080_Cyno_IC50');
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 1628.4});

%% APP siRNA Cyno IC50
variantObj = addvariant(model,'APP_siRNA_Cyno_IC50');
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 161});

%% Tominersen Mouse IC50
variantObj = addvariant(model,'Tominersen_Mouse_IC50');
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 1698.5});

%% Tominersen Cyno IC50
variantObj = addvariant(model,'Tominersen_Cyno_IC50');
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 675.02});

%% Jacifusen Mouse IC50
variantObj = addvariant(model,'Jacifusen_Mouse_IC50');
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 300.67});

%% Jacifusen Cyno IC50
variantObj = addvariant(model,'Jacifusen_Cyno_IC50');
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 1683});

%% IC50 cyno get dose
variantObj = addvariant(model,'get_dose_Cyno_IC50');
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 675.02});

%% IC50 mouse get dose
variantObj = addvariant(model,'get_dose_Mouse_IC50');
addcontent(variantObj, {'parameter', 'IC50_mRNA_Oligo_b', 'Value', 1553.4});