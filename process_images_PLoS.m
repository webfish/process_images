function process_images_PLoS
% PROCESS_IMAGES_PLOS Function for streamlined processing of images.
%   PROCESS_IMAGES_PLOS was used to prepare publication quality images
%   from raw microscopy data.
%   The function requires functioning dip_image toolbox installation.
%       <http://www.diplib.org/>
%
%   The input are four lsm files with four different fluorophore images.
%       The 405 nm excitation images contain the antibody staining.
%       The 488 nm excitation images contain the staining with anti-IgE
%           FISH probe.
%       The 594 nm excitation images contain the staining with anti-IgG
%           FISH probe.
%       The 647 nm  excitationimages contain the staining with anti-IgM
%           FISH probe.
%       
%
%   The image data is imported using tifread29jn function.
%   Spectral unmixing of the 405 nm channel from 488 nm and 594 nm channels
%       is performed. It may result in intensity dips observed in the
%       unmixed images as a result of chromatic aberrations.
%   Images are filtered by a Gaussian filter.
%   Project one or two slices containing the strongest FISH staining by
%       their averaging.
%   Threshold the 2 % of the image area with the highest intensity.
%   Select two brightest areas in the image.
%   Perform background correction for immmunoglobulin stain image and
%       suppress signal from neighboring nuclei.
%   Combine the images into a four-color merge:
%       gray background: immunoglobulin staining
%       red:             anti-IgE FISH probe
%       green:           anti-IgG FISH probe
%       blue:            anti-IgM FISH probe
%   
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details:
% <http://www.gnu.org/licenses/>.
%   
%   Copyright 2012 Jakub Nedbal, <http://www.webfish2.org>
%   $Revision: 0 $ $Date: 2012/09/27 $

%% Generate folders in not existing
if ~exist('output_images', 'dir')
    mkdir('output_images');
end
%% Generate folders in not existing
if ~exist('slices', 'dir')
    mkdir('slices');
end


%% list files and the roots to process. The full filename is the 405
%% background image and the roots following are replacing the part of the
%% first root with the subsequent roots, the slices to use and if only the
%% largest cell should be used
filelist = {'input_images/IgE-B-IgE-1.tif', '-IgE', '-488', '-594', '-647', 14, true, [0, 0; 0, 0; 0, 0], [1, 2], 0.02, [0, 118; 0, 120], 9.28; ...             % EE
            'input_images/IgE-E-IgE-1.tif', '-IgE', '-488', '-594', '-647', [12, 13], false, [0, 0; 0, 0; 0, 0], [1, 2], 0.02, [15, 127; 0, 111], 9.28; ...    % 
            'input_images/IgE-M-IgE-2.tif', '-IgE', '-488', '-594', '-647', [12, 13], false, [0, 0; 0, 0; 0, 0], [1, 2], 0.02, [5, 95; 0, 88], 9.28; ...
            'input_images/IgG1_W_IgG1-1.tif', '_IgG1', '_488', '_594', '_647', 16, false, [0, 0; 0, 0; 0, 0], [1, 2], 0.02, [8, 114; 5, 110], 9.28; ...
%            '~/biology/uscopy/2010-06-26/OB/IgG1/cutouts/IgG1_U_IgG1-2.tif', '_IgG1', '_488', '_594', '_647', 16, false, [0, 0; 2, 2; 2, 0], [1, 2], 0.005, [8, 118; 4, 112], 1.9; ...
            'input_images/IgG1-D-IgG1-1.tif', '-IgG1', '-488', '-594', '-647', [13 15], false, [0, 0; 0, 0; 0, 0], [1 2], 0.02, [3, 106; 10, 108], 1.9; ...
            'input_images/IgM-D-dapi-1.tif', '-dapi', '-488', '-594', '-647', [13], false, [0, 0; 0, 0; 0, 0], [], 0.02, [12, 112; 10, 110], 9.28};

%% overlap maximum and background maximum, image maximum
thresh = [80, 127, 255];

%% Load the files in a loop
for i = 1 : size(filelist, 1)
    for j = 2 : 5
        fname = regexprep(filelist{i, 1}, filelist{i, 2}, filelist{i, j});
        fprintf('Loading: %s\n', fname);
        %% Load the stack
        [stack, meanstack] = tiffread29jn(fname);
        XYZcoef = 8e-8;
        if j == 2
            %% Make a 4-D dip_image array to hold the entire nuclear image
            nucleus = newim(stack(1).width, stack(1).height, numel(stack), 4);
        end
        %% Convert to dip_image and crop
        for k = 0 : numel(stack) - 1 
            nucleus(:, :, k, j - 2) = dip_image(stack(k + 1).data);
        end
    end
    %% Crop the image
    nucleus = nucleus(filelist{i, 11}(1, 1) : filelist{i, 11}(1, 2), filelist{i, 11}(2, 1) : filelist{i, 11}(2, 2), :, :);

    %% Perform spectral unmixing
    nucleusU = nucleus;
    for k = filelist{i, 9}
        chX = single(nucleusU(:, :, :, k));
        chX = reshape(chX, numel(chX), 1);
        ch0 = single(nucleusU(:, :, :, 0));
        ch0 = reshape(ch0, numel(ch0), 1);
        [p, S, mu] = polyfit(ch0, chX, 1);
        subm = polyval(p, ch0, S, mu);
        subim = dip_image(reshape(subm, size(nucleusU, 2), size(nucleusU, 1), size(nucleusU, 3)));
        nucleusU(:, :, :, k) = squeeze(nucleusU(:, :, :, k)) - subim;
    end
    %% Apply gaussian filter with kernel size of one voxel to all images
    nucleusF = nucleusU;
    for k = 0 : 3
        nucleusF(:, :, :, k) = gaussf(nucleusF(:, :, :, k), 1);
    end
    %% Select desired slices and make their average projection
    nucleusS = squeeze(mean(nucleusF(:, :, filelist{i, 6}, :), [], 3));
%         %% Subtract background offset
%         nucleusB = nucleusS;
%         for k = 0 : 3
%             nucleusB(:, :, k) = backgroundoffset(squeeze(nucleusB(:, :, k)));
%         end
    %% Chromatic correction - NOT USED
    nucleusC = nucleusS;
    for k = 1 : 3
        if any(filelist{i, 8}(k, :) ~= 0)
            nucleusC(:, :, k) = shift(squeeze(nucleusC(:, :, k)), filelist{i, 8}(k, :));
        end
    end
            
    %% Threshold the image so that only 2 % of the area is shown
    nucleusT = nucleusC;
    for k = 1 : 3
        nucleusT(:, :, k) = tophat(nucleusT(:, :, k)) .* threshold(nucleusT(:, :, k), 'volume', filelist{i, 10});% .* threshold(nucleusT(:, :, 0));
    end
    %% Adjust the minimum value of the thresholded region to zero
    nucleusA = nucleusT;
    for k = 1 : 3
        mv = sort(unique(reshape(double(nucleusA(:, :, k)), prod(size(nucleusA(:, :, k))), 1)));
        mv = mv(2);         % minimum non-zero value
        nucleusA(:, :, k) = nucleusA(:, :, k) - mv;
        nucleusA(:, :, k) = nucleusA(:, :, k) - dip_image(nucleusA(:, :, k) == min(nucleusA(:, :, k)), 'double') * min(nucleusA(:, :, k));
    end
    %% Background correction for nuclear image
    %% make fake random object
    
    %nucleusA(:, :, 0) = backgroundoffset(squeeze(nucleusA(:, :, 0)));
    [~, ~, mask] = backgroundoffset(dip_image(rand(size(nucleusA, 2), size(nucleusA, 1))) * max(nucleusA(:, :, 0)) .* squeeze(threshold(nucleusA(:, :, 0))) + squeeze(~threshold(nucleusA(:, :, 0)) .*nucleusA(:, :, 0)));
    t = squeeze(nucleusA(:, :, 0));
    nucleusA(:, :, 0) = nucleusA(:, :, 0) - mean(t(mask));
    %% Find slices with sensible signal, whose std exceeds 10
    in = find(double(std(nucleusA(:, :, 1 : 3), [], [1 2])) > filelist{i, 12});
    %% Make a mask of sensible signals
    mask = newim(size(nucleusA, 1), size(nucleusA, 2), 'bin');
    for k = in(:)'
        mask = mask | squeeze(threshold(nucleusA(:, :, k), 'fixed', 0.1));
    end
    %% Generate the final immunoglobulin stain image
    Igstain = squeeze(nucleusA(:, :, 0));
    Igstain(Igstain < 0) = 0;       % get rid on negative values
    C = min([thresh(1) / max(Igstain .* mask), thresh(2) / max(Igstain)]);
                                    % scaling factor for intensity
	Igstain = C * Igstain;          % scale Igstain image
    %% Supress smaller cells if needed
    if filelist{i, 7}
        mask = threshold(Igstain);
        mask = gaussf(dip_edt(mask, [1 1], 1, 'ties'), 3);
        mask = waterseed(maxima(mask), -mask);
        %mask =  xor(mask, Igstain > 0);
        shade = bilateralf(dip_image(mask, 'sfloat'), 4) .* dip_image(Igstain > 0, 'sfloat');
        shade = max(shade) - shade;
        shade = shade - min(shade);
        shade = 1 / max(shade) * shade;
        Igstain = dip_image(xor(mask, Igstain > 0), 'sfloat') .* shade .* Igstain;
        lab = label(xor(mask, Igstain > 0) .* Igstain > 0, 1);
        meas = measure(lab, [], 'Size');
        lab(lab == find(meas.Size == max(meas.Size))) = 0;
        lab = lab == 0;
        Igstain = Igstain .* lab;
    end
    %% Generate the final FISH stain image
    FISH = nucleusA(:, :, 1 : 3);
    for k = 0 : 2
        C(k + 1) = min((thresh(3) - Igstain) ./ squeeze(FISH(:, :, k)));
    end
    for k = 0 : 2
        if any(in == k + 1)
            FISH(:, :, k) = FISH(:, :, k) * C(k + 1);
        else
            FISH(:, :, k) = FISH(:, :, k) * min(C);
        end
    end
                                    % scaling factor for intensity
    %FISH = FISH * C;
    
    % Keep the two brightest spots
    Fsum = squeeze(sum(FISH, [], 3));
    Flab = label(Fsum > 0);
    Flab = newim(size(Fsum), 'bin'); Flab(4 : end - 4, 4 : end - 4) = true;
    Flab = label(threshold(Fsum, 'volume', min(filelist{i, 10}, 0.01)) .* Flab);
    Fmeas =  measure(Flab, Fsum, 'Sum');
    [~, Fmeas] = sort(Fmeas.Sum, 'descend');
    Flab = Flab == Fmeas(1) | Flab == Fmeas(2);
    FISH = FISH .* repmat(Flab, [1, 1, size(FISH, 3)]);

    %% Make a square image
    eFISH = newim([max(size(FISH, 1), size(FISH, 2)) * [1, 1], size(FISH, 3)]);
    eFISH(round((size(eFISH, 1) - size(FISH, 1)) / 2) + (0 : size(FISH, 1) - 1), ...
          round((size(eFISH, 2) - size(FISH, 2)) / 2) + (0 : size(FISH, 2) - 1), :) = FISH;
    eIgstain = newim(max(size(Igstain, 1), size(Igstain, 2)) * [1, 1]);
    eIgstain(round((size(eIgstain, 1) - size(Igstain, 1)) / 2) + (0 : size(Igstain, 1) - 1), ...
          round((size(eIgstain, 2) - size(Igstain, 2)) / 2) + (0 : size(Igstain, 2) - 1)) = Igstain;
    eFISH = FISH;
    eIgstain = Igstain;

    %% Generate output image
    OUT{i} = eFISH + repmat(eIgstain, [1 1 3]);
    OUT{i} = resample(OUT{i}, [300 / size(OUT{i}, 1), 300 / size(OUT{i}, 2), 1]);

    %% Make a bar
    if i == 2
        OUT{i}(size(OUT{i}, 1) - 10 + round(- 2e-6 / XYZcoef * 300 / size(Igstain, 1) : 0), size(OUT{i}, 2) - 10 + (-4 : 0), :) = 255;
    end

    [~, f] = fileparts(filelist{i, 1});
    imwrite(uint8(cat(3, 255 - OUT{i}, newim(300, 300))), sprintf('output_images/cmyk%s.tif', f), 'tif', 'Compression', 'none');
    imwrite(uint8(OUT{i}), sprintf('output_images/rgb%s.tif', f), 'TIFF', 'Compression', 'none');
    noisy = nucleusS(:, :, 1 : end);
    noisy = noisy - repmat(min(noisy, [], [1 2]), [size(noisy, 1), size(noisy, 2), 1]);
    
    mask = newim(size(noisy)); mask(:, :, 0) = 255 / max(noisy);
    imwrite(uint8(noisy .* mask), sprintf('slices/R%s.tif', f), 'TIFF', 'Compression', 'none');
    mask = newim(size(noisy)); mask(:, :, 1) = 255 / max(noisy);
    imwrite(uint8(noisy .* mask), sprintf('slices/G%s.tif', f), 'TIFF', 'Compression', 'none');
    mask = newim(size(noisy)); mask(:, :, 2) = 255 / max(noisy);
    imwrite(uint8(noisy .* mask), sprintf('slices/B%s.tif', f), 'TIFF', 'Compression', 'none');
    imwrite(uint8(eIgstain / max(eIgstain) * 255), sprintf('slices/I%s.tif', f), 'TIFF', 'Compression', 'none');

    OUT{i} = joinchannels('RGB', OUT{i});

end
