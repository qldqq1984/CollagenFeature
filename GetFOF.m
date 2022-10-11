    function [ result ] = GetFOF( img1 )
    %GetFOF img1 must be two-dimensional.

    [pixelCounts, GLs] = imhist(img1);
    NM = sum(pixelCounts); % number of pixels

    % Calculate the mean
    meanGL = sum(GLs .* (pixelCounts / NM));

    varresult = 0;  % variance temp var
    skewresult = 0; % skewness temp var
    kurtresult = 0; % kurtosis temp var

    for i=0:1:length(pixelCounts)-1
      varresult = varresult + (i-meanGL)^2 * (pixelCounts(i+1)/NM);
      skewresult = skewresult + (i-meanGL)^3 * (pixelCounts(i+1)/NM);
      kurtresult = kurtresult + (i-meanGL)^4 * (pixelCounts(i+1)/NM)-3;
    end


    skewresult = skewresult * varresult ^-3; % skewness
    kurtresult = kurtresult * varresult ^-4; % kurtosis

    %energy
    energy = sum((pixelCounts / NM) .^ 2);

    %entropy 
    pI = pixelCounts / NM;
    entropy1 = -sum(pI(pI~=0) .* log2(pI(pI~=0)));

    % returns the same result as above
    %entropy2 = -sum(pI .* log2(pI+eps))
    %entropy3 = entropy(img1)

    result = [meanGL, varresult, skewresult, kurtresult, energy, entropy1];

    % calculate features using the matlab-optimized way; currently broken 
    % (different results)

    % Subtract the mean from the gray levels
    % m0 = pixelCounts - meanGL;
    % diffSquared = m0 .^ 2;
    % % variance
    % variance = sum(diffSquared .* (pixelCounts / NM))
    % 
    % %skewness
    % diffThird = m0 .^ 3;
    % skewness = variance ^ -3 * sum(diffThird .* (pixelCounts / NM))
    % 
    % %kurtosis
    %  diffFour = m0 .^ 4;
    %  kurtosis = variance ^ -4 * sum(diffFour .* ((pixelCounts / NM) - 3))
    end