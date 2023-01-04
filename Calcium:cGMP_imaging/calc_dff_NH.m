
function [dff, raw] = calc_dff_NH(Mr, Mask,n);

% n number of ROIs, 
%Mask cells for ROI defined with x and y based on ginput function


%Getting mean of lowest 100 pixel intensities/frame for bckgrdsub

% background = [];
% [x,y,z] = size(Mr);
% for j = 1:z
%     frame = Mr(:,:,j);
%     sortedpixels = sort(frame(:));
%     min100 = sortedpixels(1:100);
%     meanmin100 = mean(min100);
%     background(j) = meanmin100;
% end


for i=1:n;
    
    ind=[];raw_F=[];obgf=[];cbgf=[];fresult=[];gof=[];output=[];
    
    %so here find is used to get all the masked pixels (the only non zeros
    %in the Mask matrix)
    
    %copies the Mask pixels x times where x is # of frames such that the
    %indices are available for every frame in the rest of the analysis
    %so ind becomes a mask pixel by frame # 
ind = repmat(find(Mask{i}),[1 size(Mr,3)]);
ind = ind + repmat((0:(size(Mr,3)-1))*(size(Mr,1)*size(Mr,2)),[nnz(Mask{i}) 1]);
ind = ind(:);
raw_F = double(Mr(ind));
raw_F = reshape(raw_F,nnz(Mask{i}),size(Mr,3));
raw_F = mean(raw_F).';

obgf = raw_F;

%below is some leftover code if you want to get delta F by a fit to
%baseline (good if you have photobleaching) rather than just subtracting by
%baseline

% t is time in frames
% t = (0:(size(Mr,3)-1)).';
% [fresult,gof,output] = fit(t,raw_F,'exp2','Normalize','on');
%  [fresult,gof,output] = fit(t,raw_F(t),'exp2','Normalize','on');


% cbgf est le vecteur pour chaque frame du fit a 2 exp
% cbgf = feval(fresult,t);
% 
% % normalise a 1
% cbgf = cbgf/cbgf(1);
% tr = raw_F./cbgf;

% have user select left and right x limits within which is baseline (should
% be selecting whatever is minimum before your activity signal)
figure;plot(raw_F);
baseline = ginput;
time1(i) = round(baseline(1,1));
time2(i) = round(baseline(2,1));
clear baseline;
close all

raw(:,i)= raw_F;

% Old code
% baseline is taken from chosen timepoints

% dff(:,i) = (tr/mean(tr(time1(i):time2(i)))-1)*100;

% dff(:,i) = (bgf/mean(bgf(time1(i):time2(i)))-1)*100);

base = mean(raw_F(time1(i):time2(i)));

%the actual deltaF/F equation:
%non background sub version:
dff(:,i) = (raw_F-base)./(base)*100;

%background sub version:
%dff(:,i) = (raw_F-base)./(base-background')*100;

end;
