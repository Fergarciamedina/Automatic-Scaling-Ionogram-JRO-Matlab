function ProcessedImage = ionogramProcessing(IonogramData)
    
    % Intensity value reduction and 150km deletion
        imgDen = Denoising_1(IonogramData.totalpowerData, IonogramData.NoisetotalpowerData, IonogramData.RangeData);
    
    % Laplacian Filter
        imgLaplacian = LaplacianFilter(imgDen);
        
    % Eliminate outliners pixels
        [Projx, Projy] = Projection(imgLaplacian);
        
        var_ = AreaFiltering(imgLaplacian,Projy, 'y');
        imgOutliners = AreaFiltering(var_,Projx, 'x');  
    
    % Thresholding global 1
     %   imgThresh = ThreshImg_Global(imgOutliners, 0.45);
    
    % Morphological operation
        imgMorphOp = MorphOp(imgOutliners);
        
    ProcessedImage = imgMorphOp;
    
    % Thresholding global 2
        % imgThresh = ThreshImg_Global2(imgOutliners);

end

