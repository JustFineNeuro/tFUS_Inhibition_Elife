%% Figure 3
%% (1) the line connected and box plots for the top 2 SSDs

%n per group
N=[21, 19, 15];

%ordered: no-tFUS, Go-tFUS, Stop-tfUS
mu75G1=[0.44,0.475,0.324];
mu95G1=[0.578,0.615,0.385];

mu75G2=[0.42,0.437,0.39];
mu95G2=[0.565,0.605,0.53];

mu75G3=[0.415,0.451,0.415];
mu95G3=[0.63,0.59,0.602];



se75G1=[0.052,0.048,0.040];
se95G1=[0.037,0.047,0.043];

se75G2=[0.045,0.048,0.044];
se95G2=[0.048,0.061,0.049];

se75G3=[0.051,0.056,0.061];
se95G3=[0.054,0.081,0.063];

mksize=8;
close all

figure(1); 
subplot(131)
hold on
errorbar([.9 1.9],[mu75G1(1),mu95G1(1)],[se75G1(1),se95G1(1)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','k','Color','k','MarkerEdgeColor','k')
errorbar([1 2],[mu75G1(2),mu95G1(2)],[se75G1(2),se95G1(2)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','b','Color','k','MarkerEdgeColor','b')
errorbar([1.1,2.1],[mu75G1(3),mu95G1(3)],[se75G1(3),se95G1(3)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','g','Color','k','MarkerEdgeColor','g')
axis square
box off
title('rIFG')
xlim([0.8 2.2])
ylim([0.25 0.7]);


subplot(132)
hold on
errorbar([.9 1.9],[mu75G2(1),mu95G2(1)],[se75G2(1),se95G2(1)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','k','Color','k','MarkerEdgeColor','k')
errorbar([1 2],[mu75G2(2),mu95G2(2)],[se75G2(2),se95G2(2)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','b','Color','k','MarkerEdgeColor','b')
errorbar([1.1,2.1],[mu75G2(3),mu95G2(3)],[se75G2(3),se95G2(3)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','g','Color','k','MarkerEdgeColor','g')
axis square
box off
title('S1')
xlim([0.8 2.2])
ylim([0.25 0.7]);


subplot(133)
hold on
errorbar([.9 1.9],[mu75G3(1),mu95G3(1)],[se75G3(1),se95G3(1)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','k','Color','k','MarkerEdgeColor','k')
errorbar([1 2],[mu75G3(2),mu95G3(2)],[se75G3(2),se95G3(2)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','b','Color','k','MarkerEdgeColor','b')
errorbar([1.1,2.1],[mu75G3(3),mu95G3(3)],[se75G3(3),se95G3(3)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','g','Color','k','MarkerEdgeColor','g')
axis square
box off
title('Control')
xlim([0.8 2.2])
ylim([0.25 0.7]);



figure(1); 
subplot(131)
hold on
errorbar([.9 1.9],[mu75G1(1),mu95G1(1)],[se75G1(1),se95G1(1)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','k','Color','k','MarkerEdgeColor','k')
errorbar([1 2],[mu75G1(2),mu95G1(2)],[se75G1(2),se95G1(2)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','b','Color','k','MarkerEdgeColor','b')
errorbar([1.1,2.1],[mu75G1(3),mu95G1(3)],[se75G1(3),se95G1(3)],'Marker','o','MarkerSize',mksize,'LineWidth',1,'MarkerFaceColor','g','Color','k','MarkerEdgeColor','g')
axis square
box off
title('rIFG')
xlim([0.8 2.2])
ylim([0.25 0.7]);



%% (2) SSRT boxplots, too.
jitter_width = 0.25;
colorsdrop(1,:)=[0 0 0];
colorsdrop(2,:)=[143	38	217]/255;
colorsdrop(3,:)=[255	83	105]/255;

load('SSRTPLOT.mat')
close all
figure(2)

%% GROUP 1
for gpidx=1:3
    gp=find(SSRT_PLOT(:,end)==gpidx);

    jitter_pts=((rand(1,length(gp))-0.5)*jitter_width)';
    jitter_pts=jitter_pts;
    x=[ones(length(gp),1),2*ones(length(gp),1),3*ones(length(gp),1)]+jitter_pts;
    y=SSRT_PLOT(gp,1:3);
    y(:,2)=y(:,2)-5;


    subplot(1,3,gpidx);
    hold on
    for i=1:length(x)
        plot(x(i,:),y(i,:),'k-','LineWidth',1,'Color',[.5, .5, .5, 0.6])
    end


    for i=1:length(x)
        for j=1:3
            scatter(x(i,j),y(i,j),50,'MarkerFaceColor',colorsdrop(j,:),'MarkerEdgeColor',colorsdrop(j,:),'MarkerFaceAlpha',1);

        end
    end
    xlim([0.7 3.3])
    ylim([100 400])
    axis square
    box off
end






