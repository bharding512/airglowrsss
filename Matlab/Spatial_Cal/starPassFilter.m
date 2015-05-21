function filtered = starPassFilter(data,siz,rthresh)
% Take a series of images and return an image that accentuates features
% that move like stars.
%
% Inputs:
%     data (MxNxP) - a 3-D matrix that is constructed by stacking
%     individual images.  Specifically, the i^th image is stored in
%     data(:,:,i).
%     siz (scalar) - the size of the 3D window used in the filtering.  siz
%     is measured in pixels and should be bigger than the size of a star
%     and bigger than the number of pixels a star travels in one frame.
%     (Typical: 2).
%     rthresh (scalar) - the maximum brightness change expected in a star
%     pixel from one frame to the next.  This threshold is used to identify
%     like pixels between images.  The higher the threshold, the less noise
%     will be removed. (Typical: 200)
% Outputs:
%     filtered (MxNxP) - a 3-D matrix in the same format as "data".  Note
%     that the "borders" of this matrix are unchanged.  This means that for
%     each image, a border of size "siz" is not affected, and that the
%     first siz and last siz images are not affected.
%
% Example Usage:
%     % This code is an example of how to use starPassFitler
%     count = 55; % number of frames to use
%     tmpdata = readtif(sprintf('../Narayan/044/6300_%04i.tif',1));
%     dk = readtif('../Narayan/044/DARK_0001.tif'); % load dark frame
%     siz = size(tmpdata);
%     data = zeros(siz(1),siz(2),count); % initialize 3D matrix
%     for i = 1:count
%         % load individual images
%         fn = sprintf('../Narayan/044/6300_%04i.tif',i);
%         d = readtif(fn) - dk; % subtract dark frame
%         d(d<0) = 0; % ignore negative values
%         data(:,:,i) = d; % store in matrix
%     end
%     windowSize = 2;
%     thresh = 200;
%     filtered = starPassFilter(data,windowSize,thresh);
%
% Future work on this function:
%   1) Generalize this so there aren't any thresholds, only gradual
%      fall-offs.
%   2) Handle edge cases better (in space and time)
%   3) Instead of this ad-hoc method, use a Fourier perspective.


dsize = [siz,siz,siz]; % [i,j,time]

deli = -dsize(1):dsize(1);
delj = -dsize(2):dsize(2); 
delk = [-dsize(3):-1, 1:dsize(3)]; % don't include the current time
% Create 3-D mask that will discard the middle pixel (deli=delj=0) for all
% times (i.e. for all k)
mask = true(length(deli),length(delj),length(delk));
mask((length(deli)+1)/2,(length(delj)+1)/2,:) = 0;

% Create histogram and find the intensity threshold below which 99.9% of
% the pixels lie.  Project the image down onto this 99.9% range.
[h,b] = hist(data(:),10000);
relcs = cumsum(h)/sum(h);
thresh = b(sum(relcs<0.999));
idx = data(:)>thresh; %indices of pixels that are too bright
% now project
data(idx) = thresh;

datarmFPN = data;
[s1,s2,s3] = size(data);
for i = dsize(1)+1 : s1-dsize(1)
    for j = dsize(2)+1 : s2-dsize(2)
        for k = dsize(3)+1 : s3-dsize(3)
            ivec = i + deli;
            jvec = j + delj;
            kvec = k + delk;            
            
            windowpixels = data(ivec,jvec,kvec);
            % make sure to discard this pixel
            windowpixels = windowpixels(mask);
            
            % compute relative difference between this pixel and the
            % neighborhood
            rdist = abs(windowpixels-data(i,j,k));
            % construct similarity matrix
            similar = rdist < rthresh;
            % Collapse similarity matrix along the time dimension.
            ijwindowsum = sum(similar,3);
            % Remove a pixel if, during the interval, it correlates with
            % too few neighborhood pixels.  This avoids correlation with
            % close-by fixed-pattern noise.
            N = sum(sum(logical(ijwindowsum)));
            if N <= 1
                % replace this pixel with the median of its neighborhood
                datarmFPN(i,j,k) = median(windowpixels(:));                
            end
        end
    end
    fprintf('%i / %i\n',i-dsize(1),s1-2*dsize(1));
end

filtered = datarmFPN;


