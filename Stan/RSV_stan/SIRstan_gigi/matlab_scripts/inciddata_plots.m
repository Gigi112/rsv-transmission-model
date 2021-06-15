figure
subplot(3,2,1)
bar(cases.kolkata./population.kolkata*1e5)
title('Kolkata')
set(gca,'XTickLabel',{'<2y','2-4y','5y','10y','20y','30y','40y','50+y'})

subplot(3,2,2)
bar(cases.jakarta./population.jakarta*1e5)
title('Jakarta')
set(gca,'XTickLabel',{'<2y','2-4y','5y','10y','20y','30y','40y','60+y'})

subplot(3,2,3)
bar(cases.kibera./population.kibera*1e5)
title('Kibera')
ylabel('Incidence rate (per 100,000)')
set(gca,'XTickLabel',{'<2y','2-4y','5y','10y','18y','35y','50+y'})

subplot(3,2,4)
bar(cases.lwak./population.lwak*1e5)
title('Lwak')
set(gca,'XTickLabel',{'<2y','2-4y','5y','10y','18y','35y','50+y'})

subplot(3,2,5)
bar(cases.kalamapur./population.kalamapur*1e5)
title('Kalamapur')
set(gca,'XTickLabel',{'<2y','2-4y','5y','50+y'})

subplot(3,2,6)
bar(cases.delhi./population.delhi*1e5)
title('Delhi')
set(gca,'XTickLabel',{'<2y','2-4y','5y','13y','20y','40+y'})

%%
incid_site(:,1)=[cases.lwak(1:4); sum(cases.lwak(5:6)); cases.lwak(7)]./[population.lwak(1:4); sum(population.lwak(5:6)); population.lwak(7)]*1e5;
incid_site(:,2)=[cases.jakarta(1:4); sum(cases.jakarta(5:7)); cases.jakarta(8)]./[population.jakarta(1:4); sum(population.jakarta(5:7)); population.jakarta(8)]*1e5;
incid_site(:,3)=[cases.kolkata(1:4); sum(cases.kolkata(5:7)); cases.kolkata(8)]./[population.kolkata(1:4); sum(population.kolkata(5:7)); population.kolkata(8)]*1e5;
incid_site(:,4)=[cases.kibera(1:4); sum(cases.kibera(5:6)); cases.kibera(7)]./[population.kibera(1:4); sum(population.kibera(5:6)); population.kibera(7)]*1e5;
incid_site(:,5)=[cases.kalamapur(1:3); cases.kalamapur(3); 0; cases.kalamapur(4)]./[population.kalamapur(1:3); population.kalamapur(3:4); population.kalamapur(4)]*1e5;
incid_site(:,6)=cases.delhi./population.delhi*1e5;

figure
bar(incid_site)
legend({'Lwak','Jakarta','Kolkata','Kibera','Kalamapur','Delhi'})
set(gca,'XTickLabel',{'<2y','2-4y','5y','13y','20y','40+y'})


%%
figure
for j=1:6
subplot(3,2,j)
bar([inciddata(j).cases./inciddata(j).popn Ifit1_varR0(j).byage']*1e5)
title(inciddata(j).setting)
set(gca,'XTickLabel',inciddata(j).agecat)
end

%%
figure
for j=1:4
subplot(2,2,j)
y1=bar([inciddata2(j).cases./inciddata2(j).popn Ifit2_dr1(j).byage']*1e5);
title(inciddata2(j).setting)
set(gca,'XTickLabel',inciddata2(j).agecat)
set(y1(1),'FaceColor','b')
set(y1(2),'FaceColor','r')
if j==1
    ylabel('Incidence rate (per 100,000)')
elseif j==3
    xlabel('Age group')
end
end

%%
[kappa_est,kappa_int]=regress(log(pfit1_varR0(2,:)'./(1-pfit1_varR0(2,:)')),[lambda_varR0(1,:)' ones(6,1)])

x=.000001:.000001:.01;
figure
hold on
scatter(lambda_varR0(1,:)',pfit1_varR0(2,:)')
plot(x,1./(1+exp(-kappa_est(1)*x-kappa_est(2))),'r')
text(lambda_varR0(1,:)'+.0005,pfit1_varR0(2,:)',sites)
%plot(x,1./(1+exp(-kappa_int(1,1)*x-kappa_int(2,1))),'--r')
%plot(x,1./(1+exp(-kappa_int(1,2)*x-kappa_int(2,2))),'--r')
ylabel('Estimated reporting rate')
xlabel('Force of infection (\lambda)')

%% Maryland challenge study data
challenge1_drpar=regress(log(challenge1_p(2:end)./(1-challenge1_p(2:end))),[challenge1(2:end,1) ones(4,1)]);

figure
hold on
scatter(challenge1(:,1),challenge1_p,'ok')
plot([challenge1(:,1)'; challenge1(:,1)'],challenge1_ci','k')
plot((2.5:.1:9.5)',1./(1+exp(-challenge1_drpar(1)*(2.5:.1:9.5)'-challenge1_drpar(2))),'r','LineWidth',1)
ylabel('Probability of symptoms')
xlabel('Log10 challenge dose')
%%
figure
for j=1:4
subplot(1,4,j)
y1=bar([inciddata2(j).cases./inciddata2(j).popn Ifit2_varR0(j).byage']*1e5);
title(inciddata2(j).setting)
set(gca,'XTickLabel',inciddata(j).agecat)
set(y1(1),'FaceColor','b')
set(y1(2),'FaceColor','r')
xlim([0 length(inciddata2(j).cases)+1])
if j==1
    ylabel('Incidence (per 100,000)')
end
end
legend('Observed','Predicted')

%%
ulim=[200 100 300 50];

figure
for j=1:4
subplot(2,4,j)
hold on
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact1c(1,j).Dnovacc,2)./sum(vaccimpact1c(1,j).pop,2),'--b','LineWidth',2) %./sum(pop,2)*100000
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*(sum(vaccimpact1c(1,j).Dnovacc,2)-vaccimpact1c(1,j).popdirecteff)./sum(vaccimpact1c(1,j).pop,2),'g','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact1c(1,j).D,2)./sum(vaccimpact1c(1,j).pop,2),'r','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact1c(1,j).Dnovacc,2)./sum(vaccimpact1c(1,j).pop,2),'--b','LineWidth',2)
fill([1/52:1/52:(tspan-t0)/52 (tspan-t0)/52:-1/52:1/52 1/52]',52.18*100000*[overallimpact1c_lCI(:,j); flipud(overallimpact1c_uCI(:,j)); overallimpact1c_lCI(1,j)],[1 .8 .8],'EdgeColor',[1 .8 .8])
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*overallimpact1c_lCI(:,j),'--r','LineWidth',1)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*overallimpact1c_uCI(:,j),'--r','LineWidth',1)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*(sum(vaccimpact1c(1,j).Dnovacc,2)-vaccimpact1c(1,j).popdirecteff)./sum(vaccimpact1c(1,j).pop,2),'g','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact1c(1,j).D,2)./sum(vaccimpact1c(1,j).pop,2),'r','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact1c(1,j).Dnovacc,2)./sum(vaccimpact1c(1,j).pop,2),'--b','LineWidth',2)
%ylim([0 250]) %ulim(j)])
set(gca,'XLim',[0 (tspan-t0)/52],'XTick',5:5:(tspan-t0)/52,'XTickLabel',0:5:(tspan-t0)/52)
if j==1
ylabel({'Incidence rate'; '(per 100,000 per year)'},'FontSize',12) %,'FontSize',14)
elseif j==2
legend('No vaccination','Direct effect of vaccination','Predicted overall effect','Location','SO','Orientation','Horizontal')
end
xlabel('Year') %,'FontSize',14)
title(inciddata2(j).setting,'FontWeight','Bold','FontSize',12)

subplot(2,4,4+j)
hold on
bar([(vaccimpact1(1,j).cumincidnovacc(11,:)-vaccimpact1(1,j).cumincid(11,:))' (vaccimpact1c(1,j).cumincidnovacc(11,:)-vaccimpact1c(1,j).cumincid(11,:))'])
plot([.85:2.85; .85:2.85],[prctile(casesaverted10y(:,:,j),2.5); prctile(casesaverted10y(:,:,j),97.5)],'k')
plot([1.15:3.15; 1.15:3.15],[prctile(casesaverted10yc(:,:,j),2.5); prctile(casesaverted10yc(:,:,j),97.5)],'k')
set(gca,'XTick',1:3,'XTickLabel',{'<5y','5-20y','20+y'})
xlabel('Age group')
if j==1
ylabel({'Cases averted (per 100,000)';'over 10 years'},'FontSize',12)
legend('Routine','Routine+campaign')
end
end

%%
figure
for j=1:4
subplot(2,4,j)
hold on
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact_dr1(1,j).Dnovacc,2)./sum(vaccimpact_dr1(1,j).pop,2),'--b','LineWidth',2) %./sum(pop,2)*100000
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*(sum(vaccimpact_dr1(1,j).Dnovacc,2)-vaccimpact_dr1(1,j).popdirecteff)./sum(vaccimpact_dr1(1,j).pop,2),'g','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact_dr1(1,j).D,2)./sum(vaccimpact_dr1(1,j).pop,2),'r','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact_dr1(1,j).Dnovacc,2)./sum(vaccimpact_dr1(1,j).pop,2),'--b','LineWidth',2)
fill([1/52:1/52:(tspan-t0)/52 (tspan-t0)/52:-1/52:1/52 1/52]',52.18*100000*[overallimpact_dr1_lCI(:,j); flipud(overallimpact_dr1_uCI(:,j)); overallimpact_dr1_lCI(1,j)],[1 .8 .8],'EdgeColor',[1 .8 .8])
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*overallimpact_dr1_lCI(:,j),'--r','LineWidth',1)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*overallimpact_dr1_uCI(:,j),'--r','LineWidth',1)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*(sum(vaccimpact_dr1(1,j).Dnovacc,2)-vaccimpact_dr1(1,j).popdirecteff)./sum(vaccimpact_dr1(1,j).pop,2),'g','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact_dr1(1,j).D,2)./sum(vaccimpact_dr1(1,j).pop,2),'r','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact_dr1(1,j).Dnovacc,2)./sum(vaccimpact_dr1(1,j).pop,2),'--b','LineWidth',2)
%ylim([0 250]) %ulim(j)])
if j==1
ylabel({'Incidence rate'; '(per 100,000 per year)'},'FontSize',12) %,'FontSize',14)
elseif j==2
legend('No vaccination','Direct effect of vaccination','Predicted overall effect','Location','SO','Orientation','Horizontal')
%ylim([0 800])
end
xlabel('Year') %,'FontSize',14)
title(inciddata2(j).setting,'FontWeight','Bold','FontSize',12)

subplot(2,4,4+j)
hold on
bar(vaccimpact_dr1(1,j).cumincidnovacc(11,:)-vaccimpact_dr1(1,j).cumincid(11,:))
plot([1:3; 1:3],[prctile(casesaverted_dr1_10y(:,:,j),2.5); prctile(casesaverted_dr1_10y(:,:,j),97.5)],'k')
set(gca,'XTick',1:3,'XTickLabel',{'<5y','5-20y','20+y'})
xlabel('Age group')
if j==1
ylabel({'Cases averted (per 100,000)';'over 10 years'},'FontSize',12)
end
end


