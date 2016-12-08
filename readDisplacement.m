[dat,str] = xlsread('F:\surge_project\displacements.xlsx',1,'A2:T22')

%%
nv = logical((isnan(dat(:,12))==0))
xi = dat(nv,12)
x = dat(nv,13:15)
v = dat(nv,17:19)
x_t = cumsum([xi,x],2)
x_t_LT = x_t(:,3:4)

x_t66 = []
nv66 = logical((isnan(dat(:,4))==0))
x66 = dat(nv66,1:3)
x_t66(:,4) = dat(nv66,4)
x_t66(:,3) = x_t66(:,4)-x66(:,3)
x_t66(:,2) = x_t66(:,3)-x66(:,2)
x_t66(:,1) = x_t66(:,2)-x66(:,1)

x_t66_LT = x_t66(:,3:4)

%%
close all
figure('units','normalized','outerposition',[0 0 1 1])
% subplot(1,2,1)
cm = colormap('autumn')

for i = 1:length(x_t_LT(:,1))
    for j = 1:length(x_t_LT(1,:))
        hold on
        h{j} = plot(x_t_LT(i,j),x_t_LT(i,1),'ko','markerfacecolor',cm(j*12,:))
    end
end
grid on
ylabel('Centerline distance (km) of initial point (2013-248)')
xlabel('Centerline distance (km)')      
% legend([h{1} h{2} h{3} h{4}],{'2013-September','2014-September','2015-August','2016-September'},'location','northwest')
ylim([10 32])
xlim([12 30])


%%% last two
% subplot(1,2,2)
cm = colormap('winter')

for i = 1:length(x_t66_LT(:,1))
    for j = 1:length(x_t66_LT(1,:))
        hold on
        h{j} = plot(x_t66_LT(i,j),x_t66_LT(i,1),'ko','markerfacecolor',cm(j*12,:))
    end
end
grid on
ylabel('Centerline distance (km) of initial point (2013-248)')
xlabel('Centerline distance (km)')      
legend([h{1} h{2} ],{'1966-September','1967-August'},'location','northwest')
ylim([10 32]) 
xlim([12 30])

% %%% with monthly vel
% subplot(1,2,2)
% cm = colormap('winter')
% 
% for i = 1:length(x_t66(:,1))
%     for j = 1:length(x_t66(1,:))
%         hold on
%         h{j} = plot(x_t66(i,j),x_t66(i,1),'ko','markerfacecolor',cm(j*12,:))
%     end
% end
% grid on
% ylabel('Centerline distance (km) of initial point (2013-248)')
% xlabel('Centerline distance (km)')      
% legend([h{1} h{2} h{3} h{4}],{'1951','1966-August','1966-September','1967-August'},'location','northwest')
% ylim([10 32]) 


%% get some velocity at the midpoint

vel_66 = dat(nv66,11)
x_t66Mid = (x_t66_LT(:,1)+x_t66_LT(:,2))./2

close all
f1 = figure
h1 = plot(x_t66Mid,vel_66,'ko','markerfacecolor',[1 0.4 0],'markersize',8)
grid on


vel_15 = dat(nv,19)
x_tMid = (x_t_LT(:,1)+x_t_LT(:,2))./2
hold on
h2 = plot([6.96;x_tMid],[3.7;vel_15],'k^','markerfacecolor',[0 0 1],'markersize',8)

xlabel('Centerline midpoint distance (km)')
ylabel('Yearly average velocity during peak surge (m/d)')
legend([h1 h2],{'Sep. 1966 - Aug. 1967','Aug. 2015 - Sept. 2016'},'location','southwest')
set(gca,'fontsize',18)
savePDFfunction(f1,'F:\surge_project\code\figures\yearlyVel_comp')
















