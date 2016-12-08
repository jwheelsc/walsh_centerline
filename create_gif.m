function [] = create_gif(datType,fSpeed)
% this is the the folder where your images are located
% clear variables
% pnum = '16'
% 
% datType = ['pairs_' pnum '\accel_' pnum]
% fSpeed = 4

folderStr = ['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\strain\' datType]

imStr = cellstr(ls(folderStr))
imStr(1:2)=[]

% here is the output filname
fOut= ['C:\Users\cromp\Desktop\laptop_surge_project\code walsh\figures\strain\' datType '.gif']

% make an array of times that you want the dates the be separated by (seconds)
%%
for i = 1:length(imStr)

    v = imStr{i}
    year1(i) = str2num(v(13:16))
    day1(i) = str2num(v(18:20))
%     keyboard
end
jDay = ((year1-2013)*365)+day1
[jDay,elS] = sort(jDay)

delDay = jDay(2:end)-jDay(1:end-1)

timeArr = [delDay 32]*(1/(16*fSpeed))


    
%%
% here is a loop where I stitch the images together
for i = 1:length(imStr)
    % read in the image
    [A,cmap] = imread([folderStr '\' imStr{elS(i)}],'jpg');
    B(:,:,1) = fliplr([A(:,:,1)]');
    B(:,:,2) = fliplr([A(:,:,2)]');
    B(:,:,3) = fliplr([A(:,:,3)]');
    cmap = cmap';
%     map = map';
    % define the position of the text box
%     position = [1300 100];
    % insert the text box
%      A = insertText(A,position,imNStr{i},'FontSize',24,'BoxColor','red','BoxOpacity',0.4);
    % convert the image to RGB
    [Im cmap] = rgb2ind(B, 256);
    % append and write the png to a .gif
	if i == 1;
		imwrite(Im,cmap,fOut,'gif','LoopCount',Inf,'DelayTime',timeArr(i));
	else
		imwrite(Im,cmap,fOut,'gif','WriteMode','append','DelayTime',timeArr(i));
	end
end


%% change the file names so that the 
for i =1:length(elS)
    oldName = imStr{elS(i)}
    movefile([folderStr '\' oldName],[folderStr '\d' num2str(jDay(i)) '_' oldName])
end