%% Allometric scaling plots of oligo parameters across species
% Standalone script that plots IC50, oligonucleotide half-life and protein
% half-life versus body weight (mouse/cyno/human) on log-log axes for the
% BIIB080, Mivelsiran, Tominersen and Jacifusen programs, and exports the
% figures. Parameter values are hard-coded in the matrices below.

clear all;close all;clc;

weights = [0.025,3.8,70];
species = ["Mouse","Cyno","Human"];
bw = [0.025,3.8,70];
oligo_hl = [33.13,36.22,36.22;
NaN,118,118;    
29.12,18.52,24.16;
57.6,57.6,57.6];

% Without delay
IC50 = [117.23,1628.4,162.68;
    NaN,161,113.75;
    1698.5,675.02,826.34;
    300.67,1683,281.46];

% Delay
% IC50 = [68.42,3214,22.93;
%     NaN,148.2,22.16;
%     125,675.02,826.34;
%     300.67,1683,281.46];

% Without delay
prot_hl = [10 16 114.44;
    NaN 4.15 3.09;
    1 1 50;
    1 1 52.5];

% prot_hl = [16, 16, 104.21;
%     NaN, 2.76, 4.93;
%     1, 1, 50;
%     1, 1, 52.5];

color_scale = ["-or","-ob","-om","-ok"];

for ii=1:1

    subplot(1,1,ii)
    loglog(bw,IC50(1,:),'-.or','MarkerFaceColor','r',MarkerSize=8,LineWidth=1.5);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([1e1,1e4]);

    subplot(1,1,ii)
    loglog(bw,IC50(2,:),'--ob','MarkerFaceColor','b',MarkerSize=8,LineWidth=1.5);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([1e1,1e4]);

    subplot(1,1,ii)
    loglog(bw,IC50(3,:),':om','MarkerFaceColor','m',MarkerSize=8,LineWidth=2);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([1e1,1e4]);

    subplot(1,1,ii)
    loglog(bw,IC50(4,:),'-ok','MarkerFaceColor','k',MarkerSize=8,LineWidth=1.5);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([1e1,1e4]);

    legend(["BIIB080","Mivelsiran","Tominersen","Jacifusen"],"Location","south");


end
xticks([0.025,3.8,70]);
ylabel('IC50 (nM)')
fontsize(gcf,20,'points')
legend box off
exportgraphics(gcf,'IC50_final.png','Resolution',600);

figure
for ii=1:1

    subplot(1,1,ii)
    loglog(bw,oligo_hl(1,:),'-.or','MarkerFaceColor','r',MarkerSize=8,LineWidth=1.5);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([1e1,1e3]);

    subplot(1,1,ii)
    loglog(bw,oligo_hl(2,:),'--ob','MarkerFaceColor','b',MarkerSize=8,LineWidth=1.5);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([1e1,1e3]);

    subplot(1,1,ii)
    loglog(bw,oligo_hl(3,:),':om','MarkerFaceColor','m',MarkerSize=8,LineWidth=2);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([1e1,1e3]);

    subplot(1,1,ii)
    loglog(bw,oligo_hl(4,:),'-ok','MarkerFaceColor','k',MarkerSize=8,LineWidth=1.5);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([1e1,1e3]);

    % legend(["BIIB080","Mivelsiran","Tominersen","Jacifusen"],"Location","southeast");

    legend off

end

xticks([0.025,3.8,70]);
ylabel('Oligo half life (day)')
fontsize(gcf,20,'points')
exportgraphics(gcf,'oligo_hl_final.png','Resolution',600);

% legend box off

figure
for ii=1:1

    subplot(1,1,ii)
    loglog(bw,prot_hl(1,:),'-.or','MarkerFaceColor','r',MarkerSize=8,LineWidth=1.5);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([5e-1,5e2]);

    subplot(1,1,ii)
    loglog(bw,prot_hl(2,:),'--ob','MarkerFaceColor','b',MarkerSize=8,LineWidth=1.5);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([5e-1,5e2]);

    subplot(1,1,ii)
    loglog(bw,prot_hl(3,:),':om','MarkerFaceColor','m',MarkerSize=8,LineWidth=3);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([5e-1,5e2]);

    subplot(1,1,ii)
    loglog(bw,prot_hl(4,:),'-ok','MarkerFaceColor','k',MarkerSize=8,LineWidth=1.5);  
    hold on;
    xlim([1e-2,1e2]);
    ylim([5e-1,5e2]);

    % legend(["BIIB080","Mivelsiran","Tominersen","Jacifusen"],"Location","southeast");

    legend off

end
xticks([0.025,3.8,70]);
ylabel('Protein half life (day)')
fontsize(gcf,20,'points')
exportgraphics(gcf,'prot_hl.png','Resolution',600);

% legend box off