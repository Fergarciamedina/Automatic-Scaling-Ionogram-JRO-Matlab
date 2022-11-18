%% Plotear y guardar imagenes - O Mode Power
clear;
close all;
clc;
cd('D:\Cursos_2020_2\Tesis\Data\Codes')

IonogramaPrueba = rgb2gray(imread('imgPrueba.png'));
IonogramaPruebaInverseValues = uint8(255) - IonogramaPrueba;
IonogramaPruebaInverseValues = mat2gray(IonogramaPruebaInverseValues);

figure(1); imagesc(IonogramaPruebaInverseValues);
ylabel('Altura virtual (km)'); xlabel('Frecuencia (kHz)'); title('Ionograma de prueba'); colormap gray;
%surf(double(OmodepowerData));

ncolumnsXnrows = size(IonogramaPruebaInverseValues);
IonogramaPruebaProjMeth = IonogramaPruebaInverseValues;

%%
heightScaled = linspace(ncolumnsXnrows(1),1,1000);

% Getting rid of the traces under 150 mts and above 600 mts

for i=1:length(heightScaled)
    if (i >= 155) && (i <= 595)
 
    else
        for x = 1: ncolumnsXnrows(2)
            IonogramaPruebaProjMeth(round(heightScaled(i)), x) = 0;
        end
    end
end


figure(2); imagesc(IonogramaPruebaProjMeth); colormap gray;
ylabel('Altura virtual (km)'); xlabel('Frecuencia (kHz)'); title('Getting rid of the traces under 150 mts and above 600 mts'); colormap gray;

%% 
proyy=zeros(ncolumnsXnrows(1),1); % Para y
proyx=zeros(ncolumnsXnrows(2),1); % Para x


for y = 1:ncolumnsXnrows(1)
    for x = 1: ncolumnsXnrows(2)
        proyy(y) = proyy(y)+IonogramaPruebaProjMeth(y,x);
    end
end

for x = 1:ncolumnsXnrows(2)
    for y = 1: ncolumnsXnrows(1)
        proyx(x) = proyx(x)+IonogramaPruebaProjMeth(y,x);
    end
end
figure(4)
imagesc(IonogramaPruebaProjMeth)

figure(5)
plot(proyx);
title('Proyección X'); 

figure(6)
plot(proyy)
title('Proyección Y');



proyyFilt=zeros(ncolumnsXnrows(1),1); % Para y
proyxFilt=zeros(ncolumnsXnrows(2),1); % Para x

%% a(n) para X Este usan en el paper!!!!!
% a(n) para X
ax=proyx;
Ax=mean(ax);
k=1.7;
thresholdx = k*Ax;
OmodepowerDtFiltered1=IonogramaPruebaProjMeth;
% a(n) -> n = x 1:452 -> length(proyx) -> length(a)

for x = 1:ncolumnsXnrows(2)
    for y = 1: ncolumnsXnrows(1)
        if ax(x)<thresholdx
            OmodepowerDtFiltered1(y,x)=0;
        end
    end
end

for x = 1:ncolumnsXnrows(2)
    for y = 1: ncolumnsXnrows(1)
        proyxFilt(x) = proyxFilt(x)+OmodepowerDtFiltered1(y,x);
    end
end

figure(6)
subplot(1,2,1)
plot(proyx)
title('Proyección X');

subplot(1,2,2)
plot(proyxFilt)
title('Proyección X - Filtered');

figure(7)
subplot(1,2,1)
imagesc(IonogramaPruebaProjMeth)
title('Before Filt')
subplot(1,2,2)
imagesc(OmodepowerDtFiltered1)
title('After Filt')
%% a(n) para Y Este NOPE
% a(n) para Y
% ay=proyy;
% Ay=mean(ay);
% k=2;
% 
% thresholdy = k*Ay;
% 
% OmodepowerDtFiltered1=IonogramaPruebaProjMeth;
% 
% % a(n) -> n = x 1:452 -> length(proyx) -> length(a)
% for x = 1:ncolumnsXnrows(2)
%     for y = 1: ncolumnsXnrows(1)
%         if ay(y)<thresholdy
%             OmodepowerDtFiltered1(y,x)=0;
%         end
%     end
% end
% 
% for y = 1:ncolumnsXnrows(1)
%     for x = 1: ncolumnsXnrows(2)
%         proyyFilt(y) = proyyFilt(y)+IonogramaPruebaProjMeth(y,x);
%     end
% end
% 
% figure(8)
% subplot(1,2,1)
% plot(proyy)
% title('Proyección Y');
% 
% subplot(1,2,2)
% plot(proyyFilt)
% title('Proyección Y - Filtered');
% figure(9)
% subplot(1,2,1)
% imagesc(IonogramaPruebaProjMeth)
% title('Before Filt')
% subplot(1,2,2)
% imagesc(OmodepowerDtFiltered1)
% title('After Filt')
%% Dilation

SE = strel('rectangle',[3 8]); % lo cambie de 3x3 a 4x6 para que se parezca al paper
IonogramDil = imdilate(OmodepowerDtFiltered1,SE);

figure, imagesc(OmodepowerDtFiltered1)
title('Ionograma')
saveas(gcf,'IonogramaPaperPrueba.png')

figure, imagesc(IonogramDil), title('Ionograma dilatado')
saveas(gcf,'DilatacionIonogramaPaperPrueba.png')

%% connected component labeling

IonogramCCL = bwlabel(IonogramDil); % Connected Componet Labeling
IonogramCCLFiltered = IonogramCCL;

rCenter=0;
% cCenter=0;
 
for ii = min(min(IonogramCCL)):max(max(IonogramCCL))
    
    [r, c] = find(IonogramCCL==ii);
    
    rCenter = round((min(r)+max(r))/2);
    % cCenter = round((min(c)+max(c))/2);

    if (rCenter > round(heightScaled(100))) && (rCenter < round(heightScaled(550)))
    else
        IonogramCCLFiltered(r,c)=0; % Delete isolated
    end
end

% rCenterScaled = linspace(ncolumnsXnrows(1),1,1000);

figure(9); 
subplot(1,2,1); 
figure, imagesc(IonogramCCL); title('Ionograma - CCL')
saveas(gcf,'CCLIonogramaPaperPrueba.png')
subplot(1,2,2);imagesc(IonogramCCLFiltered); title('Ionograma - CCL Filtered')

%%
IonogramCCL = bwlabel(IonogramDil); % Connected Componet Labeling

IonogramCCLFiltered = IonogramCCL;

% ----

[r, c] = find(IonogramCCL==2);

IonCenter = IonogramCCL(((min(r)):max(r)) , (min(c):max(c)));
figure; imagesc(IonCenter)

% floor((size(IonogramCCLFiltered(r,c))+1)/2) % center

IonogramCCLFiltered(r,c)=0; % Delete isolated

figure(10); 
subplot(1,2,1); imagesc(IonogramCCL); title('Ionograma - CCL')
subplot(1,2,2); imagesc(IonogramCCLFiltered); title('Ionograma - CCL Filtered Manual')

%% Overlap
IonogramCCLFiltered=mat2gray(IonogramCCLFiltered);
IonDone = OmodepowerDtFiltered1 - IonogramCCLFiltered;

for x = 1:ncolumnsXnrows(2)
    for y = 1: ncolumnsXnrows(1)
        if IonDone(y,x) <0.50
            IonDone(y,x)=0;
        end
    end
end

% zzz = imoverlay(IonogramCCLFiltered,OmodepowerDtFiltered1);

figure(11); 
subplot(1,2,1); imagesc(IonogramaPruebaInverseValues); title('Ionograma'); colormap('gray');
subplot(1,2,2); imagesc(IonDone); title('Ionograma - Final'); colormap('gray');
%% 

% plot(proyx); axis('off');
% saveas(gcf,'1.png')
imshow(imread('2.png'))

zzz1 = rgb2gray(imread('1.png'));
zzz2 = rgb2gray(imread('2.png'));

% Get size of existing image A.
[rowsA, colsA, numberOfColorChannelsA] = size(zzz2);
% Get size of existing image B.
[rowsB, colsB, numberOfColorChannelsB] = size(zzz1);
% See if lateral sizes match.
if rowsB ~= rowsA || colsA ~= colsB
    % Size of B does not match A, so resize B to match A's size.
    zzz1 = imresize(zzz1, [rowsA colsA]);
end


% imshow(zzz2)
err = immse(zzz2, zzz1);
figure; imshow(zzz2)
figure; imshow(zzz1)
%%

plot(zzz1(:,1))
