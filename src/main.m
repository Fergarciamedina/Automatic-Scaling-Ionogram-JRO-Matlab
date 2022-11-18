clc; clear; clf, close all;

Folders.srcPath  = pwd;

Folders.Functions = [Folders.srcPath, '\Functions\'];
Folders.NGI = uigetdir(pwd); 

fileList = dir(fullfile(Folders.NGI, '/*.ngi'));

Folders.imgProc = [Folders.NGI, '\Imagenes_Procesadas\'];

if ~exist(Folders.imgProc, 'dir')
   mkdir(Folders.imgProc)
end

for i = 1:length(fileList)

    % Get data
    cd(Folders.Functions)
    IonogramData = GettingRelevantData(fileList(i).name, Folders);
    
    % Processing
    ProcessedImage  = ionogramProcessing(IonogramData);
    
    CustomPlot(IonogramData, ProcessedImage, Folders.imgProc, fileList(i).name)
    
    cd(Folders.srcPath)
end
clear ProcessedImage
%%
figure('visible','on');
var_ = GettingRelevantData(fileList(282).name, Folders);
imagesc(var_.totalpowerData)

