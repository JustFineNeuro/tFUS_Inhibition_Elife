%% P100 errorbar plot
%Same as in JASP
load('p100.mat')
SS=mean(p100(:,[1 3 5]),2);
US=mean(p100(:,[2 4 6]),2);

z=1.96; %95% confidence interval

CI_SS_pos=mean(SS)+z*(std(SS)./sqrt(length(SS)));
CI_SS_neg=mean(SS)-z*(std(SS)./sqrt(length(SS)));

CI_US_pos=mean(US)+z*(std(US)./sqrt(length(US)));
CI_US_neg=mean(US)-z*(std(US)./sqrt(length(US)));

close all
figure; 
hold on
for i=1:length(SS)
plot(1:2,[SS(i),US(i)],'color', [.5 .5 .5], 'linewidth', 1.5)
end
scatter(ones(length(SS),1),SS,70,'Marker','o','MarkerFaceColor',[230 97 0]/255,'MarkerEdgeColor',[230 97 0]/255);
scatter(2*ones(length(US),1),US,70,'Marker','o','MarkerFaceColor',[93 58 155]/255,'MarkerEdgeColor',[93 58 155]/255);

errorbar([.96],mean([SS]),[CI_SS_neg],[CI_SS_pos],'Color','k','LineWidth',3)
errorbar([2.04],mean([US]),[CI_US_neg],[CI_US_neg],'Color','k','LineWidth',3)
scatter(.96,mean([SS]),100,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')
scatter(2.04,mean([US]),100,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')

xlim([0.88 2.15])
axis square