%%
clc
clear variables
close all
vpA = []
sos = 27
%%%this is 33 for lucania



pairNum = 16
pnum = num2str(pairNum)

%%% read in the landsat 8 image
% iL8 = imread('F:\surge_project\LC80630172015133LGN00_B8.TIF');
load('cSeas.mat')
load('walshContC.mat')
% load('keep_i_6364.mat')
seasB = length(seas(:,1))
load('WalshCenterDots.mat')
load('keep_i_walsh.mat')
load('WalshMaskDots.mat')
centerline = maskDots;

cC = centerlineCont;
% cC = centerLHod;
diffD = cC(2:end,:)-cC(1:end-1,:);
centerD = [0;sqrt(sum(diffD.^2,2))];
ccD = cumsum(centerD);

folderStr_63 = 'C:\Users\cromp\Desktop\laptop_surge_project\p063_r017\'
lsNames_63 = cellstr(ls([folderStr_63 '*_0' pnum '_*.nc']));
lsNamesTif_63 = cellstr(ls([folderStr_63 '*_0' pnum '_*.tif']));
for i = 1:length(lsNames_63)

    fName_63{i} = [folderStr_63 lsNames_63{i}];
    fNameTif_63{i}  = [folderStr_63 lsNamesTif_63{i}]; 
    
end
folderStr_62 = 'C:\Users\cromp\Desktop\laptop_surge_project\p062_r017\'
lsNames_62 = cellstr(ls([folderStr_62 '*_0' pnum '_*.nc']));
lsNamesTif_62 = cellstr(ls([folderStr_62 '*_0' pnum '_*.tif']));
for i = 1:length(lsNames_62)

    fName_62{i} = [folderStr_62 lsNames_62{i}];
    fNameTif_62{i}  = [folderStr_62 lsNamesTif_62{i}]; 
    
end


% 
% folderStr_64 = 'D:\Documents\Courses\surge_project\velocity_dec05\p064_r017\'
% lsNames_64 = cellstr(ls([folderStr_64 '*_0' pnum '_*.nc']));
% lsNamesTif_64 = cellstr(ls([folderStr_64 '*_0' pnum '_*.tif']));
% for i = 1:length(lsNames_64)
% 
%     fName_64{i} = [folderStr_64 lsNames_64{i}];
%     fNameTif_64{i}  = [folderStr_64 lsNamesTif_64{i}]; 
%     
% end
% 
% fName = [fName_62,fName_63,fName_64];
% fNameTif = [fNameTif_62,fNameTif_63,fNameTif_64];

% fName = [fName_63];
% fNameTif = [fNameTif_63];

fName = [fName_63,fName_62];
fNameTif = [fNameTif_63,fNameTif_62];

for i = 1:length(fName)

    varN = fName{i};
    namArr{i} = varN([50:57]+sos)
    
end
[blank,elSrt] = sort(namArr);
fNameA = fName(elSrt);
fNameTifA = fNameTif(elSrt);

loopArr = 1:length(fName)
% loopArr(not_i)=[] 

%%
close all
f1 = figure('units','normalized','outerposition',[0 0 1 1])
fs = 14
count = 1

for ii = 1:length(keep_i)
    i = keep_i(ii)
% for i= 1:length(loopArr) 
% for i = 62
    close all
    f1 = figure
    fName = fNameA{i};
    fNameTif = fNameTifA{i}; 

    vv = ncread(fName,'vv_masked');
    vv = vv';
%     ncdisp(fName,'vv_masked')
%% read in the log10 geotif
    [X,cmap] = geotiffread(fNameTif);

    xArr = linspace(cmap.XWorldLimits(1),cmap.XWorldLimits(2),length(X(1,:,1)));
    yArr = linspace(cmap.YWorldLimits(2),cmap.YWorldLimits(1),length(X(:,1,1)));
     

%all the walsh
    xRange = [495000 560000];
    yRange = [6740000 6770000];

    
    elXLw = find(abs(xArr-xRange(1))<1500);
    xl = elXLw(end);
    elXUp = find(abs(xArr-xRange(2))<1000);
    xu = elXUp(end);
    elYLw = find(abs(yArr-yRange(1))<1000);
    yl = elYLw(end);
    elYUp = find(abs(yArr-yRange(2))<1000);
    yu = elYUp(end);
%% read in the band 8 geotiff
%     [X8,cmap] = geotiffread('E:\surge_project\LC80620172013248LGN00_B8.TIF');
%     %%
%     xArr8 = linspace(cmap.XWorldLimits(1),cmap.XWorldLimits(2),length(X8(1,:)));
%     yArr8 = linspace(cmap.YWorldLimits(2),cmap.YWorldLimits(1),length(X8(:,1)));
% 
%     xRange = [530000 557000];
%     yRange = [6770000 6799000];
% 
%     elXLw = find(abs(xArr8-xRange(1))<1000);
%     xl8 = elXLw(end);
%     elXUp = find(abs(xArr8-xRange(2))<1000);
%     xu8 = elXUp(end);
%     elYLw = find(abs(yArr8-yRange(1))<1000);
%     yl8 = elYLw(end);
%     elYUp = find(abs(yArr8-yRange(2))<1000);
%     yu8 = elYUp(end);
%     
%     subplot(1,2,1)
%     imagesc(xArr8,yArr8,X8)
%     xlim([xArr8(xl8) xArr8(xu8)])
%     ylim([yArr8(yl8) yArr8(yu8)])
%     colormap gray
%     set(gca,'xdir','reverse')

%     
    subplot(1,2,1)
    imagesc(xArr,yArr,(vv))
    xlim([xArr(xl) xArr(xu)])
    ylim([yArr(yl) yArr(yu)])
    c = colorbar;
    c.Label.String = 'm/d';
    caxis([0 12]);
%     alpha(ih,'0.1')
    ylabel('UTM northing')
    xlabel('UTM easting')
    set(gca,'xdir','reverse')
    set(gca,'fontsize',fs)
%     hold on 
%     plot(OL2(:,1),OL2(:,2),'r')
%     hold on
%     plot(centerline(:,1),centerline(:,2),'r+')
% % %     


    clML = [];
    for j = 1:length(centerline(:,1))

        px = centerline(j,1);
        py = centerline(j,2);
        p1 = find(abs(xArr-px)<150.2);
        p1 = p1(end);
        clML(j,1) = p1;
        p2 = find(abs(yArr-py)<150.2);
        clML(j,2) = p2;

    end
%     clMLu = unique(clML,'rows');
%     hold on
    clMLu = clML;
    hold on
    
    for j = 1:length(clMLu(:,1))
        xEl = clMLu(j,1);
        yEl = clMLu(j,2);
        xC = xArr(xEl);
        yC = yArr(yEl);
        subplot(1,2,1)
        hold on
        plot(xC,yC,'r+')
        vp(j) = vv(yEl,xEl);
        
        dcC = cC-repmat([xC,yC],length(cC(:,1)),1);
        dC = sqrt(sum(dcC.^2,2));
        min_dCel = find(dC==min(dC));
        dvp(j) = ccD(min_dCel(1));
       
    end
    sp1.NextPlot = 'ReplaceChildren';
    sp2 = subplot(1,2,2)
    day1 = str2num(fName([55:57]+sos));
    fracD = day1/365;
    elSB = ceil(fracD*seasB);
    elSBA(count) = elSB;
    colorB = seas(elSB,end:-1:1);
    plot(dvp/1000,vp,'ko','markerfacecolor',colorB)
    ylim([0 17 ])
    xlim([0 70])
    grid on
    title(['path ' fName([38:40]+sos) ' row ' fName([42:44]+sos) ' -- from ' fName([50:53]+sos) '-' fName([55:57]+sos) ' to ' fName([59:62]+sos) '-' fName([64:66]+sos) '---i = ' num2str(i)])
    sp2.NextPlot = 'ReplaceChildren';

    set(gca,'fontsize',fs)
% 
    vpA(count,:) = vp; 
    dvpS(count,:) = dvp;
%     pause(.1)
    ylabel('Speed (m/d)')
    xlabel('Centerline distance (km)')
    cm2 = colormap(sp2,seas(:,end:-1:1));
    cb2 = colorbar;
%     cb2.TickLabels = []
    cb2.Ticks = [0 0.25 0.5 0.75 1]
    cb2.TickLabels = {'January','April','July','October','January'}
%     set(ax, 'xdir','reverse')
% 
    savePDFfunction(f1,['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\vel_mask\' fName([38:66]+sos) ])

    count = count+1


end

return
%%
namArr = namArr(keep_i)

for i = 1:length(namArr)
    
    v = namArr{i}
    year1(i) = str2num(v(1:4))
    day1(i) = str2num(v(6:8))

end

jDay = ((year1-2013)*365)+day1
[jDay,elS] = sort(jDay)



[x1,els] = sort(dvpS(1,:))
dvpSs = dvpS(:,els)
vpAs = vpA(:,els)
save(['vel_walsh_mask_' pnum '.mat'],'dvpSs','vpAs','namArr')

%%
close all
f2 = figure
ih = imagesc(x1,jDay,vpAs)

xls = get(gca,'xlim')
hold on
plot([xls(1) xls(2)],[365 365],'r')
hold on
plot([xls(1) xls(2)],[365*2 365*2],'r')
hold on
plot([xls(1) xls(2)],[365*3 365*3],'r')
t = text(0.1e4,120,'2013')
t.Color = 'r'
t = text(0.1e4,365+30,'2014')
t.Color = 'r'
t = text(0.1e4,365*2+30,'2015')
t.Color = 'r'
t = text(0.1e4,365*3+30,'2016')
t.Color = 'r'
ylabel('Day since Jan 01 - 2013')
xlabel('Centerline distance (m)')
cb = colorbar
cb.Label.String = 'Speed (m/d)'
caxis([0 12])
set(gca,'fontsize',18)  
title(['Walsh glacier velocity from ' pnum ' day pairs'])
savePDFfunction(f2,['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\vel_spaceTime_mask'])
% 
%%

close all
f3 = surf(jDay,x1,vpAs')
xlabel('Day')
ylabel('Centerline distance (m)')
zlabel('Speed (m/day)')
zlim([0 5])
grid on
view([-170 45])
hold on
yl = get(gca,'ylim')
plot([365 365],yl,'r')
hold on
plot([365*2 365*2],yl,'r')
hold on
plot([365*3 365*3],yl,'r')
savefig(['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\vel_spaceTime'])
%% just velocity 
close all

delDay = jDay(2:end)-jDay(1:end-1)

timeArr = [delDay pairNum]*(1/(16*4))
f3=figure
for i = 1:length(keep_i)
    ii = keep_i(i)
    fName = fNameA{ii};
    colorB = seas(elSBA(i),end:-1:1);
    plot(dvpSs(i,:),vpAs(i,:),'ko','markerfacecolor',colorB)
    ylim([0 17])
    xlim([0 70000])
    grid on
%     colormap 'winter'
    colormap(seas(:,end:-1:1))
    cb2 = colorbar
    cb2.Ticks = [0 0.25 0.5 0.75 1]
    cb2.TickLabels = {'January','April','July','October','January'}
    xlabel('Centerline distance (m)')
    ylabel('Speed (m/d)')
    set(gca,'fontsize',18)
    title(['path ' fName([38:40]+sos) ' row ' fName([42:44]+sos) ' -- from ' fName([50:53]+sos) '-' fName([55:57]+sos) ' to ' fName([59:62]+sos) '-' fName([64:66]+sos)])
    savePDFfunction(f3,['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\vel_Walsh_dt\' fName([38:66]+sos)])
%     pause(timeArr(i))

    
end


%% just velocity continuous

close all
delDay = jDay(2:end)-jDay(1:end-1)

timeArr = [delDay pairNum]*(1/(16*4))
f3=figure
for i = 1:length(keep_i)
    ii = keep_i(i)
    fName = fNameA{ii};
    hold on
    colorB = seas(elSBA(i),end:-1:1);
    plot(dvpSs(i,:),vpAs(i,:),'ko','markerfacecolor',colorB)
    ylim([0 17])
    xlim([0 70000])
    grid on
%     colormap 'winter'
    colormap(seas(:,end:-1:1))
    cb2 = colorbar
    cb2.Ticks = [0 0.25 0.5 0.75 1]
    cb2.TickLabels = {'January','April','July','October','January'}
    xlabel('Centerline distance (m)')
    ylabel('Speed (m/d)')
    set(gca,'fontsize',18)
    title(['path ' fName([38:40]+sos) ' row ' fName([42:44]+sos) ' -- from ' fName([50:53]+sos) '-' fName([55:57]+sos) ' to ' fName([59:62]+sos) '-' fName([64:66]+sos)])
    savePDFfunction(f3,['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\vel_Walsh_dtCont\' fName([38:66]+sos)])
%     pause(timeArr(i))
    
end



%% velocity with compression
delDay = jDay(2:end)-jDay(1:end-1)

timeArr = [delDay pairNum]*(1/(16*4))
f4=figure()
for i = 1:length(keep_i)
    ii = keep_i(i)
    fName = fNameA{ii};
    colorB = seas(elSBA(i),end:-1:1);
    
    yyaxis left
    hold off
    x = dvpSs(i,:)
    y1 = vpAs(i,:)
    y = y1(isnan(y1)==0)
    x = x(isnan(y1)==0)
    [x,el]=sort(x)
    y = y(el)
    plot(x,y,'ko','markerfacecolor',colorB)
    ylim([0 15])
    xlim([0 70000])
    grid on
    colormap(seas(:,end:-1:1))
    cb2 = colorbar
    cb2.Ticks = [0 0.25 0.5 0.75 1]
    cb2.TickLabels = {'January','April','July','October','January'}
    xlabel('Centerline distance (m)')
    ylabel('Speeed (m/d)')
    set(gca,'fontsize',18)
    title(['path ' fName([38:40]+sos) ' row ' fName([42:44]+sos) ' -- from ' fName([50:53]+sos) '-' fName([55:57]+sos) ' to ' fName([59:62]+sos) '-' fName([64:66]+sos)])
%     savePDFfunction(f3,['F:\surge_project\code\figures\vel_32_dt\' fName(38:66)])
    

    yyaxis right 
    hold off
    dvds = -(diff(y)./diff(x))
    xmid = x(1:end-1)+(diff(x)/2)
    plot(xmid,dvds,'ko-','markerfacecolor','r')
    ylim([-2e-3  3e-3])
    hold on
    plot(get(gca,'xlim'),[0 0],'k')
    yar.NextPlot = 'ReplaceChildren';
    savePDFfunction(f4,['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\vel_compress\' fName([38:66]+sos)])
end

%% going to try to get acceleration

close all
delDay = jDay(2:end)-jDay(1:end-1)

timeArr = [delDay pairNum]*(1/(16*4))
f3=figure
for i = 1:length(keep_i)-1
    ii = keep_i(i)
    fName = fNameA{ii};
    colorB = seas(elSBA(i),end:-1:1);
    y = (vpAs(i+1,:) - vpAs(i,:))./delDay(i)
    plot(dvpSs(i,:),y,'ko','markerfacecolor',colorB)
    xlim([0 70000])
    grid on
%     colormap 'winter'
    colormap(seas(:,end:-1:1))
    cb2 = colorbar
    cb2.Ticks = [0 0.25 0.5 0.75 1]
    cb2.TickLabels = {'January','April','July','October','January'}
    xlabel('Centerline distance (m)')
    ylabel('Acceleration (m/d^{2})')
    set(gca,'fontsize',18)
    ylim([-2e-1 2e-1])
    title(['path ' fName([38:40]+sos) ' row ' fName([42:44]+sos) ' -- from ' fName([50:53]+sos) '-' fName([55:57]+sos) ' to ' fName([59:62]+sos) '-' fName([64:66]+sos)])
    savePDFfunction(f3,['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\accel\' fName([38:66]+sos)])
%     pause(timeArr(i))
    pause(0.5)
    
end

%% going to try to get acceleration and vel in the same subplot



delDay = jDay(2:end)-jDay(1:end-1)

timeArr = [delDay pairNum]*(1/(16*4))
f3=figure
for i = 1:length(keep_i)-1
    ii = keep_i(i)
    fName = fNameA{ii};
    
    subplot(2,1,1)
    colorB = seas(elSBA(i+1),end:-1:1);
    plot(dvpSs(i+1,:),vpAs(i+1,:),'ko','markerfacecolor',colorB)
    ylim([0 17])
    xlim([0 70000])
    grid on
%     colormap 'winter'
    colormap(seas(:,end:-1:1))
    cb2 = colorbar
    cb2.Ticks = [0 0.25 0.5 0.75 1]
    cb2.TickLabels = {'January','April','July','October','January'}
%     xlabel('Centerline distance (m)')
    ylabel('Speed (m/d)')
    set(gca,'fontsize',18)
    title(['path ' fName([38:40]+sos) ' row ' fName([42:44]+sos) ' -- from ' fName([50:53]+sos) '-' fName([55:57]+sos) ' to ' fName([59:62]+sos) '-' fName([64:66]+sos)])

    subplot(2,1,2)
    colorB = seas(elSBA(i),end:-1:1);
    y = (vpAs(i+1,:) - vpAs(i,:))./delDay(i)
    plot(dvpSs(i,:),y,'ko','markerfacecolor',colorB)
    xlim([0 70000])
    grid on
%     colormap 'winter'
    colormap(seas(:,end:-1:1))
    cb2 = colorbar
    cb2.Ticks = [0 0.25 0.5 0.75 1]
    cb2.TickLabels = {'January','April','July','October','January'}
    xlabel('Centerline distance (m)')
    ylabel('Acceleration (m/d^{2})')
    set(gca,'fontsize',18)
    ylim([-2e-1 2e-1])
%     title(['path ' fName(38:40)  ' row ' fName(42:44) ' -- from ' fName(50:53) '-' fName(55:57) ' to ' fName(59:62) '-' fName(64:66)])
    savePDFfunction(f3,['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\velANDaccel\' fName([38:66]+sos)])
%     keyboard

end