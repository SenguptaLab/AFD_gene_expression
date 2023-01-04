function M = multitiff2M(tiff,framesel);
%framesel should be the set of frames you want to grab (i.e. 1:361)
tiffinfo = imfinfo(tiff);
% infos sur image, struct dont la taille est le nb de frames

numframes = length(tiffinfo);
%nb de frames

M = zeros(tiffinfo(1).Height,tiffinfo(1).Width,['uint' num2str(tiffinfo(1).BitDepth)]);

for frame=framesel;
    curframe = im2uint16(imread(tiff,frame));
    M(:,:,frame-framesel(1)+1) = curframe;
   
end
