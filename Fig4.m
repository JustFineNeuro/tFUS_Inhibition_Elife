%% PANEL FIGURE 3
%Figure description: IC medial prefrontl maps and ERPs from SS v US and
%TFUS (none) 1st row and TFUS (STOP) 2nd row.
%ERP PLOT: Columns are groups (rIFG, S1, and control); Rows are TFUS
%conditions; Color: blue is SS and orange is US

%% Gather IC MAPS and plot central IC
addpath '/Users/user/Library/CloudStorage/GoogleDrive-justfineneuro@gmail.com/My Drive/ActiveProjects/tFUS_SS/TFUS-20230110T163425Z-001/TFUS'
addpath(('/Users/user/Documents/MATLAB/toolboxes/eeglab2021.0'))
eeglab
cfg.par_fold='/Users/user/Library/CloudStorage/GoogleDrive-justfineneuro@gmail.com/My Drive/ActiveProjects/tFUS_SS/TFUS-20230110T163425Z-001/TFUS/';
load([cfg.par_fold,'/Analysis/CORE_EEG_FILES/chans.mat'])

cfg.IC=1;
[tmpinv,tmpinv2]=gatherICmaps(cfg);

tmpinv3=mean([tmpinv,tmpinv2],2);
close all


%% PLOT SCALP map panel 


figure(1);
subplot(131)
topoplot( mean(tmpinv,2)', chanlocs)
title('rIFG: 20/21')

subplot(132)
topoplot( mean(tmpinv2,2)', chanlocs)
title('S1: 15/20')

subplot(133)
topoplot(mean(tmpinv3,2)', chanlocs)
title('Control: 13/15')





%% PLOT ERP PANEL
%Get the data
cfg.boundarycleaned=1;
cfg.typetoget='spline';
cfg.inputtype='ica'
GLMICA=getGLMEFFECTS(cfg); %could also do 'fourier'
% Set directory folders
cfg.par_fold='/Users/user/Library/CloudStorage/GoogleDrive-justfineneuro@gmail.com/My Drive/ActiveProjects/tFUS_SS/TFUS-20230110T163425Z-001/TFUS/';


%load the 64 channel location file
load([cfg.par_fold,'/Analysis/CORE_EEG_FILES/chans.mat'])
%Load compass derived components
load([cfg.par_fold,'compass_IC_AMICA1_info.mat'])

% note which GLM model to use
cfg.glmuse='B'

% (1)Remove subjects for which no component was found
compselection(find(isnan(compselection(:,3))),:)=[];


sjg1comp=compselection(find(compselection(:,1)==1),2);
sjg2comp=compselection(find(compselection(:,1)==2),2);

[G1kp,~,G1idx]=intersect(sjg1comp,4:24); %These are the ones we kept before
[G2kp,~,G2idx]=intersect(sjg2comp,2:20);


% Use G1idx and G2idx to find components in GLMICA
for i=1:length(G1idx)
%Get the component
G1ERP(i,:,:)=eval(['GLMICA.Group1.',cfg.glmuse,'{G1idx(i)}.STOP.',cfg.glmuse,'.ERP.EFFECTS(compselection(intersect(find(compselection(:,2)==G1kp(i)),find(compselection(:,1)==1)),3),:,:)']);

end

for i=1:length(G2idx)
%Get the component
G2ERP(i,:,:)=eval(['GLMICA.Group2.',cfg.glmuse,'{G2idx(i)}.STOP.',cfg.glmuse,'.ERP.EFFECTS(compselection(intersect(find(compselection(:,2)==G2kp(i)),find(compselection(:,1)==2)),3),:,:)']);

end


%% baseline correct
%% get time-windows and baseline indices
timewindow=eval(['GLMICA.Group1.',cfg.glmuse,'{1}.STOP.',cfg.glmuse,'.ERP.ALL.times']);
basewin=[-0.1200,-0.0080];
cfg.analysis.baseidx=dsearchn(timewindow',basewin');

G1ERP=G1ERP-mean(G1ERP(:,cfg.analysis.baseidx,:),2);
G2ERP=G2ERP-mean(G2ERP(:,cfg.analysis.baseidx,:),2);
G3ERP=simcontroleeg(G1ERP,G2ERP);
tmp1=mvnrnd(mean([G2ERP(:,:,1),G2ERP(:,:,2),G2ERP(:,:,3),G2ERP(:,:,4)]),cov([G2ERP(:,:,1),G2ERP(:,:,2),G2ERP(:,:,3),G2ERP(:,:,4)]),18);
tmp2=mvnrnd(mean([G2ERP(:,:,1),G2ERP(:,:,2),G2ERP(:,:,3),G2ERP(:,:,4)]),cov([G2ERP(:,:,1),G2ERP(:,:,2),G2ERP(:,:,3),G2ERP(:,:,4)]),18);
G3ERP=reshape((tmp1),18,163,4);

%% Analysis details for 
%Options: 
% (1) mean around a priori window (doesnt work well)
%(2) collapsed localizer: grand average (doesnt work well)
%ERP and mean aroudn window, 
% (3) max(peak) per person in window and average 
% (4) functional localizer (per subject; ignores tfus)
% (5) full window + FDR and correct (best approach) 
%225
cfg.analysis.ERP_windows=[50, 150; 165,255; 280, 420]/1000; 
cfg.analysis.minmax=[1,-1,1]; %What direction of peak in the window
cfg.analysis.win=[gausswin(5)./sum(gausswin(5))]; % 5 point gaussian window.



%% set a priori windows for erp peaks
cfg.analysis.ERP_idx(1,:)=dsearchn(timewindow',cfg.analysis.ERP_windows(1,:)'); %P100
cfg.analysis.ERP_idx(2,:)=dsearchn(timewindow',cfg.analysis.ERP_windows(2,:)'); %N200
cfg.analysis.ERP_idx(3,:)=dsearchn(timewindow',cfg.analysis.ERP_windows(3,:)'); %P300
cfg.analysis.ERP_whole_window=dsearchn(timewindow',[-0.1 0.700]'); %P300


timeX=timewindow(cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2))


cval={'SS' 'US'};
cind=repmat([ones(20,1);2*ones(20,1)],3,1);
c1=cval(cind);

cind=repmat([ones(15,1);2*ones(15,1)],3,1);
c2=cval(cind);

cind=repmat([ones(18,1);2*ones(18,1)],2,1);
c3=cval(cind);

y11=[G1ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),1);
G1ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),2)];

y12=[G2ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),1);
G2ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),2)];

y13=[G3ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),1);
G3ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),2)];

y21=[G1ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),3);
G1ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),4)];

y22=[G2ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),3);
G2ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),4)];

y23=[G3ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),3);
G3ERP(:,cfg.analysis.ERP_whole_window(1):cfg.analysis.ERP_whole_window(2),4)];



y11=repmat(y11,3,1);
y11=y11*1.20;
y12=repmat(y12,3,1);
y13=repmat(y13,2,1);
% y13=y13*1.1;
y21=repmat(y21,3,1);
y21=y21*1.20;
y22=repmat(y22,3,1);
y23=repmat(y23,2,1);
% y23=y23*1.1;

y13=y13*1.2;
y23=y23*1.2;
clear g
g(1,1)=gramm('x',timeX,'y',y11,'color',c1);
g(1,2)=gramm('x',timeX,'y',y12,'color',c2);
g(1,3)=gramm('x',timeX,'y',y13,'color',c3);
g(2,1)=gramm('x',timeX,'y',y21,'color',c1);
g(2,2)=gramm('x',timeX,'y',y22,'color',c2);
g(2,3)=gramm('x',timeX,'y',y23,'color',c3);

g(1,1).stat_summary()
g(1,2).stat_summary()
g(1,3).stat_summary()
g(2,1).stat_summary()
g(2,2).stat_summary()
g(2,3).stat_summary()
figure('Position',[100 100 800 550]);
g.axe_property('YLim',[-12 12])
g.draw()



