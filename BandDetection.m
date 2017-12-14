clear all
clc
close all
%% parameters
pkProminence=0.012; %prominencia de intesidad para detección de bandas
allowedError=0.02; %error permitido en la construcción de árboles
file='images/gel25.jpg';
%% read
I = imread(file);
if(size(I,3)~=1)
    I=I(:,:,1:3);
end
% crop
[I2, rect] = imcrop(I); 
close all;

if(size(I2,3)~=1)
    I2p=rgb2gray(I2);
else
    I2p=I2;
end
% set saturation contrast
I3 = imadjust(I2p);

% invert image
choice = questdlg('What color are the bands?', ...
	'Black/White', ...
	'White','Black','Black');
switch choice
    case 'Black'
       I4=I3;
    case 'White'
        I4 = imcomplement(I3);
end
 
%Select pozos
fig4 = figure(4);
imshow(I4);
pause(0.00001);
frame_h = get(handle(fig4),'JavaFrame');
set(frame_h,'Maximized',1); 

title ('Draw rectangles over each lane','fontsize',18);

%Select rectangle with mouse 
%out = [xmin ymin width height]
rect(1,:) = getrect(fig4);
hold on
rectangle('Position',rect,'EdgeColor','m','LineWidth',2);
hold off

boo = 1;
n = 1;

while boo == 1
choice = questdlg('Ready?', ...
	'Drawing rectangles...', ...
	'Selection complete','Not yet','Not yet');
switch choice
    case 'Not yet'
        n = n+1;
        rect(n,:) = getrect(fig4);
        hold on
        rectangle('Position',rect(n,:),'EdgeColor','m','LineWidth',2);
        hold off
    case 'Selection complete'
%         msgbox('oki')
        boo = 0;
end
end

%% Plots selected wells
close all
m = ceil(n/6.1);%filas plot max columnas (6.1)
nn = ceil(n/m)*2;%columnas plot
maxYpk=0;

%pp = counter crop well + intensity plot
%p odd pp vales
%q even pp values
fig8=figure(8)
for pp = 1:2*n
    p = floor(pp/2)+1;
    q = pp/2;
    
if mod(pp,2) == 1
subplot(m,nn,pp)
    %x min = distance form the edge of the image to first side of the
    %rectangle
    x1 = ceil(rect(p,1));
    %x max = distance form the edge of the image to second side of the
    %rectangle (min + width)
    x2 = ceil(rect(p,1)+rect(p,3));
    %Image: From the cropped image just draw the selected well
    Ip = I4(:,x1:x2);
    [rows,columns,numColorChannels] = size(Ip);
    %conserves the selected heigth with a width of 50 pixels
    B = imresize(Ip, [rows 50]);
    imshow(B);

elseif mod(pp,2)== 0
subplot(m,nn,pp)
    x1 = ceil(rect(q,1));
    x2 = ceil(rect(q,1)+rect(q,3));
    Ip = I4(:,x1:x2);
%Convert form rgb to double the intensity image values    
ydouble = im2double(Ip(:,end-floor(size(Ip,2)*0.5)));
ydouble = -ydouble;
%Filter the noise form the image. Savitzky Golay Filter
ys = sgolayfilt(ydouble,1,17);
% ys=ydouble;
x = 1:1:length(ys);
[amp,y_pks]=findpeaks(ys,x,'MinPeakProminence',pkProminence);
if max(y_pks)>maxYpk
    maxYpk=max(y_pks);
end

end
end
table=zeros(maxYpk,n);
for pp = 1:2*n
    p = floor(pp/2)+1;
    q = pp/2;
    
if mod(pp,2) == 1
subplot(m,nn,pp)
    %x min = distance form the edge of the image to first side of the
    %rectangle
    x1 = ceil(rect(p,1));
    %x max = distance form the edge of the image to second side of the
    %rectangle (min + width)
    x2 = ceil(rect(p,1)+rect(p,3));
    %Image: From the cropped image just draw the selected well
    Ip = I4(:,x1:x2);
    [rows,columns,numColorChannels] = size(Ip);
    %conserves the selected heigth with a width of 50 pixels
    B = imresize(Ip, [rows 50]);
    imshow(B);

elseif mod(pp,2)== 0
subplot(m,nn,pp)
    x1 = ceil(rect(q,1));
    x2 = ceil(rect(q,1)+rect(q,3));
    Ip = I4(:,x1:x2);
%Convert form rgb to double the intensity image values    
ydouble = im2double(Ip(:,end-floor(size(Ip,2)*0.5)));
ydouble = -ydouble;
%Filter the noise form the image. Savitzky Golay Filter
ys = sgolayfilt(ydouble,1,17);
x = 1:1:length(ys);
%Finds the local maxima of peaks with a prominance > 0.02
[amp,y_pks]=findpeaks(ys,x,'MinPeakProminence',pkProminence);

%Plot intensity vs Distance in y
plot(ys,-x);
axis([-1 0 min(-x) max(-x)])
hold on
plot(amp,-y_pks,'ro');
hold off

%%Matrix of intensity in each well with 1%error bands
errorP=ceil(maxYpk*allowedError/2);
for index= 1:1:size(y_pks')
    ypk=y_pks(index);
    table(ypk-0,q)=1;
    for e=1:1:errorP
        table(ypk-e,q)=1;
        table(ypk+e,q)=1;
    end
end

end
end
pause(0.00001);
frame_h = get(handle(fig8),'JavaFrame');
set(frame_h,'Maximized',1); 
% set threshold

% Each colum is a different profile
matrix_profile=table;

% m number of rows n number of columns
[m,n] = size(matrix_profile);
dice_sim = zeros(n,n);
numerator=0;

%Matrix similarity by Dice coefficient
%Dice(Si,Sj)=[%]
%2 âˆ— number of common bands in lanes i and j/(number bands in lane i +
%number bands in lane j)


for j = 1:n

    for k = 1:n
     
       vv3= matrix_profile(:,j)+matrix_profile(:,k);
       numerator=length(find(vv3==2));
             
       dice_sim(j,k)=2*numerator/sum(vv3);
        

    end

end

% Distance Matrix
dist_mat = 1-dice_sim;


%Names lanes
name = 'Lane';
for i= 1:n
   names(i,:)=[name '_' num2str(i,'%02d')];  % "Lane" with num lane
end

names = cellstr(names); %cells with char data


% treeN = seqlinkage(dist_mat);
% phytreeviewer(treeN)
% phytreewrite('E:\USUARIOS/USER/Documents/Alejo/MBC/Algoritmos/nwktrees/newtree.nwk',treeN)

% %Nearest distance (single linkage method)
% phylotree_NS = seqlinkage(dist_mat,'single',names);
% view(phylotree_NS) 
% 
% %Furthest distance (complete linkage method)
% phylotree_FD = seqlinkage(dist_mat,'complete',names);
% view(phylotree_FD) 

% % Unweighted Pair Group Method Average (UPGMA, group average).
% phylotree_UPGMA = seqlinkage(dist_mat,'average',names);
% view(phylotree_UPGMA) 

% %Weighted Pair Group Method Average (WPGMA)
% phylotree_WPGMA = seqlinkage(dist_mat,'weighted',names);
% view(phylotree_WPGMA) 
% 
% %Unweighted Pair Group Method Centroid (UPGMC)
% phylotree_UPGMC = seqlinkage(dist_mat,'centroid',names);
% view(phylotree_UPGMC) 
% 
% %Weighted Pair Group Method Centroid (WPGMC)
% phylotree_WPGMC = seqlinkage(dist_mat,'median',names);
% view(phylotree_WPGMC) 

%Neigbhor Joining
phylotree_NJ = seqneighjoin(dist_mat,'equivar',names);
view(phylotree_NJ)
phytreewrite('E:\USUARIOS/USER/Documents/Alejo/MBC/Algoritmos/nwktrees/newtree.nwk',phylotree_NJ)






