%% Plotear y guardar imagenes - O Mode Power
clear;
close all;
clc;

zippedDataFrame = 'JM91J_2020276055804.ngi';

RangeData = ncread(zippedDataFrame,'Range'); % km
FreqData = ncread(zippedDataFrame,'Frequency');  %kHz
O_mode_powerData = ncread(zippedDataFrame,'O-mode_power');
NoiseO_mode_powerData = ncread(zippedDataFrame,'O-mode_noise')';
% X_mode_powerData = ncread(zippedDataFrame,'X-mode_power');
% total_powerData = ncread(zippedDataFrame,'total_power');
% total_powerData = ncread(zippedDataFrame,'total_noise');
% ncdisp(zippedDataFrame)

freqStart = ncread(zippedDataFrame, 'freq_start'); 
freqEnd = ncread(zippedDataFrame, 'freq_end');

freq = double(FreqData);
range = double(RangeData);

% Plot O_mode_powerData
figure(10);
titlelabel = strrep(zippedDataFrame,'.ngi','');
O_mode_powerDt=mat2gray(double(O_mode_powerData));
h = pcolor(freq,range,O_mode_powerDt);
set(h, 'edgecolor','none');         %get rid of edgelines
caxis([0 1])                   % specify colorbar limits
colorbar;
% gray image
colormap('gray');

set(gca,'YDir','normal') 
ylabel('Altura virtual (km)');
xlabel('Frecuencia (kHz)');
title(['O mode - Ionograma ', strrep(titlelabel,'_','\_')]);
% saveas(gcf,[titlelabel,'_O_mode_power.png'])
%figure(11)
% surf(double(O_mode_powerData));
%% Denoising

for i = 1:452
    for ii = 1:512
        zzz(ii, i) = double(O_mode_powerData(ii, i)) - double(NoiseO_mode_powerData(i));
        if zzz(ii, i) < 0 
            zzz(ii, i) = 0; 
        end
    end
end

figure(2);
titlelabel = strrep(zippedDataFrame,'.ngi','');
test=mat2gray(double(zzz));
test = medfilt2(test);
h = pcolor(freq,range,test);
set(h, 'edgecolor','none');         %get rid of edgelines
caxis([0 1])                   % specify colorbar limits
colorbar;
% gray image
colormap('gray');

set(gca,'YDir','normal') 
ylabel('Altura virtual (km)');
xlabel('Frecuencia (kHz)');
title(['O mode - Ionograma ', strrep(titlelabel,'_','\_')]);

%% Plot O_mode_powerData
total_powerData=mat2gray(double(total_powerData));
figure(12);
imagesc(freq, range, double(total_powerData), [0,1]);
%colormap('gray');
colorbar;
set(gca,'YDir','normal') 
ylabel('Altura virtual (km)');
xlabel('Frecuencia (kHz)');
title(['Total power - Ionograma ', titlelabel]);

%% nel
clear all;
close all;
clc;

zippedDataFrame = 'JM91J_2020276055804.ngi';
%ncdisp(zippedDataFrame);

RangeData = ncread(zippedDataFrame,'Range'); % km
FreqData = ncread(zippedDataFrame,'Frequency');  %kHz
timexd = ncread(zippedDataFrame,'Time');  %kHz
O_mode_powerData = ncread(zippedDataFrame,'O-mode_power');
O_mode_noise_powerData = ncread(zippedDataFrame,'O-mode_noise');
X_mode_powerData = ncread(zippedDataFrame,'X-mode_power');
total_powerData = ncread(zippedDataFrame,'total_power');

%RangeData(50) - RangeData(49)
%1RangeData(21) - RangeData(20)

%FreqData(50) - FreqData(49)
%FreqData(21) - FreqData(20)

%whos RangeData;
%whos FreqData;
% whos O_mode_powerData;
% whos X_mode_powerData;
% whos total_powerData;

%surf(double(O_mode_powerData));
%title('O_mode_powerData');
freq = double(FreqData);
range = double(RangeData);


% Plot O_mode_powerData
figure(10);
% imagesc(freq, range, double(O_mode_powerData), [0,50]);
% set(gca,'YDir','normal') 
% colorbar;
% ylabel('Altura virtual (km)');
% xlabel('Frecuencia (kHz)');
titlelabel = strrep(zippedDataFrame,'.ngi','');
% title(['O mode - Ionograma ', strrep(titlelabel,'_','\_')]);
% saveas(gcf,[titlelabel,'_O_mode_power.png'])
h = pcolor(freq,range,double(O_mode_powerData));
set(h, 'edgecolor','none');         %get rid of edgelines
%caxis([0 1])                   % specify colorbar limits
colorbar;
%colormap('gray');
%saveas(gcf,'Ionogram.png')
set(gca,'YDir','normal') 
ylabel('Altura virtual (km)');
xlabel('Frecuencia (kHz)');
title(['Total power - Ionograma ', titlelabel]);
saveas(gcf,[titlelabel,'_O_mode_power.png'])

%%
% Plot X_mode_powerData
figure(11);
imagesc(freq, range, double(X_mode_powerData), [0,50]);
set(gca,'YDir','normal') 
colorbar;
ylabel('Altura virtual (km)');
xlabel('Frecuencia (kHz)');
title(['X mode - Ionograma ', titlelabel]);

% Plot O_mode_powerData
total_powerData=mat2gray(double(total_powerData));
figure(12);
imagesc(freq, range, double(total_powerData), [0,1]);
%colormap('gray');
colorbar;
set(gca,'YDir','normal') 
ylabel('Altura virtual (km)');
xlabel('Frecuencia (kHz)');
title(['Total power - Ionograma ', titlelabel]);

%total_powerData=mat2gray(double(total_powerData));
figure(12);
h = pcolor(freq,range,total_powerData);
set(h, 'edgecolor','none');         %get rid of edgelines
caxis([0 1])                   % specify colorbar limits
colorbar;
colormap('gray');
%saveas(gcf,'Ionogram.png')
set(gca,'YDir','normal') 
ylabel('Altura virtual (km)');
xlabel('Frecuencia (kHz)');
title(['Total power - Ionograma ', titlelabel]);

mat2gray(O_mode_powerData);
max(mat2gray(O_mode_powerData));
%%
I= rgb2gray(imread('Barbara-original-image.png'));
%I = O_mode_powerDt;
F=fft2(I);


figure; imagesc(I)
figure; image(I)
figure; imagesc(I); colormap gray;
figure; imagesc(I); colormap('gray');axis('image');
figure; image(I); colormap('gray');axis('image');axis('image','off');
min(I,[],2)
min(I,[],1)
figure; imagesc(I,[11 255]); colormap('gray');axis('image');axis('image','off');
figure; imagesc(I,[50 255]); colormap('gray');axis('image');axis('image','off');
figure; imagesc(log10(abs(F))); colormap('gray');axis('image');
figure; imagesc(fftshift(log10(abs(F))+1)); colormap('gray');axis('image');