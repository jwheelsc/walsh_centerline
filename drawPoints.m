%% this section is for making a mask
% hold on
% maskDots = []
% save('WalshMaskDots.mat','maskDots')
for l = 1:10
    c = ginput(15)
    maskDots = [maskDots;c]
    plot(maskDots(:,1),maskDots(:,2),'ro')
    save('WalshMaskDots.mat','maskDots','-append')
    keyboard
end

%%

load('WalshCenterDots.mat')
hold on
plot(centerDots(:,1),centerDots(:,2),'ro')

load('walshContC.mat')
plot(centerlineCont(:,1),centerlineCont(:,2),'r')

%%
centerlineCont =[]
save('walshContC.mat','centerlineCont')


for l = 1:10
    h = imfreehand()
    centerlineCont = [centerlineCont;h.getPosition]
    plot(centerlineCont(:,1),centerlineCont(:,2),'r')
    save('walshContC.mat','centerlineCont','-append')
    keyboard
end


%%
hold on
% centerDots = []
% save('WalshCenterDots.mat','centerDots')
for l = 1:10
    c = ginput(3)
    centerDots = [centerDots;c]
    plot(centerDots(:,1),centerDots(:,2),'ro')
    save('WalshCenterDots.mat','centerDots','-append')
    keyboard
end
