%Uses MAGATanalyzer package (https://github.com/samuellab/MAGATAnalyzer),
%which also includes the LabView package for tracking the worms. Much of
%the code here is taken directly from 'starter_commands_worm.html'

%you need to have already, while in the directory where you did the LabView
%tracking, created a timing file (.tim), perhaps with
%createTimingInForFileSets.m

%load tracks (select the .bin file)
eset = ExperimentSet.fromFiles('minpts', 480);

%*Cleaning up tracks

% To stitch together tracks in an experiment (to clean up missing frames)
doc Experiment.stitchTracks

maxFrameInterval = 4; % stitch together tracks if first ended 4 or fewer frames before second started
maxDist = 8; % stitch together tracks if first ended within 8 pixels of second's start
eset.executeExperimentFunction('stitchTracks', maxFrameInterval, maxDist);
% each time the program finds tracks to stitch, it will show you a graph (in
% figure 1) of the two tracks, with the starting track in green, the
% track to be merged onto the end in red, and any other nearby tracks in gray
% in figure 2, it will play a movie of the end of track 1 followed by the
% start of track 2
% you will get a prompt asking stitch tracks y/n
% to stitch the tracks, type 'y' and hit enter;  to leave them separate 'n'
% and enter.  to see the movie again, just hit enter

%%% CLEANING UP BAD TRACKS
% To clean up short tracks, tracks where the worm just sits there, etc.

% for information on the ESetCleaner class
doc ESetCleaner

ecl = ESetCleaner();
%ecl.minSpeed = xx; %minimum average speed of the track to be kept
%ecl.minDist = xx; %minimum distance must travel from the starting point to be kept
%ecl.minPts = xx; %minimum number of points to be kept

ecl.minPts = 750;
ecl.minSpeed = 0.4;

% generate graphs showing the distributions of speed, dist, and pts and
% which would be cut
ecl.getReport (eset);

% actually do the trimming (will show you report and ask if you want to
% proceed)
ecl.clean(eset);

%load text file with assay parameters (select the *header.txt)
%Need to select comma, 15

text_filename = uigetfile('*.*');
text_file = uiimport(text_filename);
text_file = struct2cell(text_file);

%get assay dimensions
left_side = cell2mat(text_file{1}(11));
left_side = str2double(left_side(5:length(left_side)));
right_side = cell2mat(text_file{1}(13));
right_side = str2double(right_side(5:length(right_side)))+left_side;
top = cell2mat(text_file{1}(12));
top = str2double(top(5:length(top)));
bottom = cell2mat(text_file{1}(14));
bottom = str2double(bottom(5:length(bottom)))+top;
assay_width = right_side - left_side;

%goal of this loop is to get all "real" end points
end_points = [];
start_points = [];
for i = 1:length(eset.expt.track)
    i
    end_frame = [];
    
        %this part gets the starting points
        if eset.expt.track(i).startFrame <= 100
            start_position = eset.expt.track(i).pt(1).loc;
            start_points = [start_points, start_position];
        end
        
    %if track ends near last frame, get end position
        if eset.expt.track(i).endFrame > 7000
            end_frame = eset.expt.track(i).npts;
            end_position = eset.expt.track(i).pt(end_frame).loc;
            end_points = [end_points, end_position];
        end
    
        %if track ends before last frames AND ends near far side of screen, get
        %end position
        if eset.expt.track(i).endFrame <= 7000
            end_frame = eset.expt.track(i).npts;
            end_position = eset.expt.track(i).pt(end_frame).loc;
            if end_position(1) >= (right_side-100)
                end_points = [end_points, end_position];

            elseif end_position(1) <= (left_side+100)
                end_points = [end_points, end_position];

            elseif end_position(2) <= (top+100)
                end_points = [end_points, end_position];

            elseif end_position(2) >= (bottom-100)
                end_points = [end_points, end_position];
            end
        end
    
end

%get all the points and find the average starting position
end_points_x = end_points(1,:);
end_points_y = end_points(2,:);
start_points_x = start_points(1,:);
start_line_x = mean(start_points_x);
end_points_x_from_start = end_points_x - start_line_x;

%make pretty figure
figure;hold on;eset.expt.track.plotPath('displacement','r-');%plot(end_points_x_from_start,end_points_y, 'bo'); 
plot (0,0, 'g.', 'MarkerSize', 20);
xlim([-(assay_width/2),assay_width/2]); %ylim([top-20,bottom+20]);
set(gcf,'Renderer','painters');
saveas(gcf,'end_points_super','fig');

%mean absolute end position
mean_end_points_x = mean(end_points_x);

%distance from start (probably more relevant for statistics)
end_point_dist_x = end_points_x - start_line_x;
mean_end_point_dist_x = mean(end_point_dist_x);

%see all distances from start
figure;hold on; plot(end_point_dist_x,1:length(end_point_dist_x),'o');xlim([(0-(start_line_x-left_side)),(0+(right_side-start_line_x))]);

%normalize end point distances to assay width and calculate mean
end_distances_x_norm = end_point_dist_x/assay_width;
mean_end_distances_x_norm = mean(end_distances_x_norm);

%get mean speed for all tracks and mean of means
for i=1:length(eset.expt.track)
    track_speeds(i) = mean(eset.expt.track(i).dq.speed);
end
mean_track_speeds = mean(track_speeds);

filename = 'analysis_end_points';
save(filename);

