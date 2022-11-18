%% Plotear y guardar imagenes - O Mode Power
clear;
close all;
clc;

myFolder = pwd;
folderI = 'D:\Cursos_2020_2\Tesis\Data\Imagenes_ionogramas\Facil';

cd(folderI);
zippedDataFrame = 'JM91J_2020276055804.ngi';

RangeData = ncread(zippedDataFrame,'Range'); % km
FreqData = ncread(zippedDataFrame,'Frequency');  %kHz
% OmodepowerData = ncread(zippedDataFrame,'O-mode_power');
% XmodepowerData = ncread(zippedDataFrame,'X-modepower');
totalpowerData = ncread(zippedDataFrame,'totalpower');

freq = double(FreqData);
range = double(RangeData);

% Plot totalpowerData
figure(1);
titlelabel = strrep(zippedDataFrame,'.ngi','');
totalpowerDt=mat2gray(double(totalpowerData));

h = pcolor(freq,range,totalpowerDt);
set(h, 'edgecolor','none');         %get rid of edgelines
caxis([0 1])                   % specify colorbar limits
colorbar;
% gray image
%colormap('gray');

set(gca,'YDir','normal') 
ylabel('Altura virtual (km)');
xlabel('Frecuencia (kHz)');
title(['O mode - Ionograma ', strrep(titlelabel,'','\')]);
saveas(gcf,[titlelabel,'totalpower.png'])
%figure(11)
%surf(double(totalpowerData));
cd(myFolder);
%%
%totalpowerDtFlip = flip(totalpowerDt);
totalpowerDtFlipFiltered = totalpowerDt;

for y = 1:512
    for x = 1: 452
        if totalpowerDt(y,x)<=0.5
            totalpowerDtFlipFiltered(y,x) = 0;
        end
    end
end

figure(2)
imagesc(totalpowerDtFlipFiltered)
%set(h, 'edgecolor','none');         %get rid of edgelines
caxis([0 1])                   % specify colorbar limits
colorbar;
% gray image
colormap('gray');

%% 
%totalpowerDtFlip = flip(totalpowerDt);
% totalpowerDtFlip = totalpowerDtFlipFiltered;
ncolumnsXnrowstotalpowerArray = size(totalpowerDt);
totalpowerDtprefilter = totalpowerDt;
% Getting rid of the first 50 mts

for y = 1:50
    for x = 1: 452
        totalpowerDtprefilter(y,x)=0;
    end
end

proyy=zeros(ncolumnsXnrowstotalpowerArray(1),1); %512
proyx=zeros(ncolumnsXnrowstotalpowerArray(2),1); %452
for y = 1:512
    for x = 1: 452
        proyy(y) = proyy(y)+totalpowerDtprefilter(y,x);
    end
end

for x = 1:452
    for y = 1: 512
        proyx(x) = proyx(x)+totalpowerDtprefilter(y,x);
    end
end
figure(3)
imagesc(totalpowerDtprefilter)

figure(4)
plot(proyx)
title('Proyección X');
figure(5)
plot(proyy)
title('Proyección Y');
%% a(n) para X en freq
% a(n) para X
ax=proyx;
Ax=mean(ax);
k=1.05;
thresholdx = k*Ax;
totalpowerDtFiltered1=totalpowerDtprefilter;
% a(n) -> n = x 1:452 -> length(proyx) -> length(a)
for x = 1:452
    for y = 1: 512
        if ax(x)<thresholdx
            totalpowerDtFiltered1(y,x)=0;
        end
    end
end
figure(6)
imagesc(totalpowerDtFiltered1)
%% a(n) para Y
% a(n) para Y
ay=proyy;
Ay=mean(ay);
k=1.02;
% k=1.1;
% k=1.2;
thresholdy = k*Ay;

totalpowerDtFiltered1=totalpowerDtprefilter;

% a(n) -> n = x 1:452 -> length(proyx) -> length(a)
for x = 1:452
    for y = 1: 512
        if ay(y)<thresholdy
            totalpowerDtFiltered1(y,x)=0;
        end
    end
end
figure(7)
imagesc(totalpowerDtFiltered1)
% figure(8)
% imagesc(totalpowerDt)