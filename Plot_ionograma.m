%% Plotear y guardar imagenes
clear;
close all;
clc;

zippedDataFrame = 'JM91J_2020276200803.ngi';

RangeData = ncread(zippedDataFrame,'Range'); % km
FreqData = ncread(zippedDataFrame,'Frequency');  %kHz
O_mode_powerData = ncread(zippedDataFrame,'O-mode_power');
% X_mode_powerData = ncread(zippedDataFrame,'X-mode_power');
% total_powerData = ncread(zippedDataFrame,'total_power');

freq = double(FreqData);
range = double(RangeData);

% Plot O_mode_powerData
figure(10);
titlelabel = strrep(zippedDataFrame,'.ngi','');
h = pcolor(freq,range,double(O_mode_powerData));
set(h, 'edgecolor','none');         %get rid of edgelines
caxis([0 50])                   % specify colorbar limits
colorbar;
%colormap('gray');
set(gca,'YDir','normal') 
ylabel('Altura virtual (km)');
xlabel('Frecuencia (kHz)');
title(['O mode - Ionograma ', strrep(titlelabel,'_','\_')]);

% saveas(gcf,[titlelabel,'_O_mode_power.png'])
