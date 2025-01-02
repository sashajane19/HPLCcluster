%Sasha Kramer
%sasha.kramer@lifesci.ucsb.edu
%UCSB IGPMS

%%%Code to cluster HPLC pigments and plot dendrograms
%Map the directory where you will load your data:
cd /Users/skramer/Documents/UCSB/Research/Data/HPLC_Aph_Rrs/

%Load your samples (formatted here as a .mat file):
load Global_HPLC_all.mat  

%The order of columns in this matrix (called "Global_HPLC") is:
%Tchla,Tchlb,Tchlc,ABcaro,ButFuco,HexFuco,Allo,Diadino,Diato,Fuco,Perid,Zea,MVChla,DVchla,Chllide,MVChlb,DVchlb,Chlc12,Chlc3,Lut,Neo,Viola,Phytin,Phide,Pras

%Initial quality control - set all values below detection (based on NASA GSFC limits) to 0
%0.001 = Chlc3, Chlc12, Chllide, Viola, Diadino, Diato, Allo, Zea, Lut, ABcaro
indPH = find(Global_RHPLC(:,24) ~= 0 & Global_RHPLC(:,24) <= 0.003); %Phide
Global_RHPLC(indPH,24) = 0; clear indPH

indPE = find(Global_RHPLC(:,11) ~= 0 & Global_RHPLC(:,11) <= 0.003); %Perid
Global_RHPLC(indPE,11) = 0; clear indPE

indBF = find(Global_RHPLC(:,5) ~= 0 & Global_RHPLC(:,5) <= 0.002); %ButFuco
Global_RHPLC(indBF,5) = 0; clear indBF

indFU = find(Global_RHPLC(:,10) ~= 0 & Global_RHPLC(:,10) <= 0.002); %Fuco
Global_RHPLC(indFU,10) = 0; clear indFU

indNE = find(Global_RHPLC(:,21) ~= 0 & Global_RHPLC(:,21) <= 0.002); %Neo
Global_RHPLC(indNE,21) = 0; clear indNE

indPR = find(Global_RHPLC(:,25) ~= 0 & Global_RHPLC(:,25) <= 0.002); %Pras
Global_RHPLC(indPR,25) = 0; clear indPR

indHF = find(Global_RHPLC(:,6) ~= 0 & Global_RHPLC(:,6) <= 0.002); %HexFuco
Global_RHPLC(indHF,6) = 0; clear indHF

indCB = find(Global_RHPLC(:,16) ~= 0 & Global_RHPLC(:,16) <= 0.003); %MVchlb
Global_RHPLC(indCB,16) = 0; clear indCB

indDC = find(Global_RHPLC(:,14) ~= 0 & Global_RHPLC(:,14) <= 0.002); %DVchla
Global_RHPLC(indDC,14) = 0; clear indDC

indPY = find(Global_RHPLC(:,23) ~= 0 & Global_RHPLC(:,23) <= 0.003); %Phytin
Global_RHPLC(indPY,23) = 0; clear indPY

%Check percent of pigments below detection (if >80% of values are below detection, I
%remove the pigment from the cluster analysis):
for i = 1:25
    belowD = find(Global_RHPLC(:,i) <= 0.001);
    j(i) = length(belowD);
    percent(i) = 100*(j(i)./145); %need to replace the denominator with your # of samples (here, 145)
end
clear belowD i j

%First remove degradation products: Chllide,Phytin,Phide
deg = [15,23,24];
Rpigcluster1 = Global_RHPLC;
Rpigcluster1(:,deg) = [];

%Remove redundant pigs (Tchlb, Tchlc, ABcaro, MVchla):
%In my dataset, Lut, Pras, and DVchlb were also below detection >80% of the
%time so I remove them here
deg2 = [2,3,4,8,9,13,16,19,22]; %remember that the order of pigments changed when you removed degradation pigments above
Rpigcluster2 = Rpigcluster1;
Rpigcluster2(:,deg2) = [];
label2 = {'TChla','ButFuco','HexFuco','Allo','Fuco','Perid','Zea','DVchla','MVchlb','Chlc12','Chlc3','Neo','Viola'};

%Cluster pigments:
D2 = pdist(Rpigcluster2','correlation'); %correlation is the method - you can vary that here
Z2 = linkage(D2,'ward'); %ward is the method - you can vary that here

%Plot dendrogram of absolute pigment values:
figure(2),clf
h = dendrogram(Z2,'Labels',label2);
set(gca,'XTickLabelRotation',90,'fontsize',18)
set(h,'color','k','linewidth',2)
ylabel('Linkage Distance')
ax = gca;
ax.YGrid = 'on';
box on
clear ax h

%Normalize to chlorophyll-a and re-cluster:
normchl = Rpigcluster2(:,2:end)./Global_RHPLC(:,1);
normlabel = label2(2:end);

D3 = pdist(normchl','correlation'); %correlation
Z3 = linkage(D3,'ward');

%Plot dendrogram of normalized pigment values:
figure(3),clf
h = dendrogram(Z3,'Labels',normlabel);
set(gca,'XTickLabelRotation',90,'fontsize',20)
set(h,'color','k','linewidth',2)
ylabel('Linkage Distance')
ax = gca;
ax.YGrid = 'on';
box on
clear ax h D2 D3 deg deg2 labels Z2 Z3

%Calculate the cophenetic distances and a corrcoef/pvalue for the dendrogram
%Example using Spearman's rank coeff
[c,D] = cophenet(Z3,D3);
[rho,pval] = corr(D3',D','type','spearman');

%For further analysis: cluster samples - change "maxclust" based on your dendrogram results
%Need to also look at pigment concentrations in each cluster to check taxonomic meaning
C = clusterdata(normchl,'distance','correlation','linkage','ward','maxclust',4);
