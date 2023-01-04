function [T_star_temp, total_fret_amp] = T_star_loop(Image_Time,raw_temperature,dff);

%NH: Putting in some comments to help me understand HBell's code

steep_points = [];             % This is the important part for finding the first 
time_steep = [];               % sustained increase in the del_F signal.
increasing = [];

del_F = smooth(dff,20,'loess');
diff_del_F = diff(dff);

    for j = 1:length(del_F)-20
        if del_F(j+1)>del_F(j)
            increasing(end+1) = 1;
        elseif del_F(j+1)<del_F(j)
            increasing(end+1) = 0;
        end
        if max(del_F(j:j+10))<2,
            increasing(j) = nan;
        end                        % finds points in del_F where nothing rises above 2 for 10
                                   % consective time points and changes 
                                   % just the first point to nan, not all ten
                                   % points.
    end
    %output of this loop is the variable increasing
    
    inc_points = zeros(1, length(increasing)); 
    for k = 20:length(increasing)-10
        if increasing(k:k+7) == 1  % finds points in del_F where there is continuous
            inc_points(k) = 1;     % increase over at least 8 seconds. The first point
        end                        % will be the T_set_time. The lag will factor in later
        if mean(diff_del_F(k+2:k+7))<0.3  
            inc_points(k) = 0;     %This sets the slope threshold. It checks that the slope of the 
                                   %points 3-8 past the point of interest
                                   %average at least 0.3. This is a bit
                                   %conservative.
        end                        
    end
    %the output of this loop is the variable inc_points, which indicates
    %which points in increasing satisfy continous condition
    
    %just take the correct points out to a new variable
    response_times = find(inc_points==1);
    
    %I guess this is saying take the first increase if it exists, otherwise
    %wait for user input?
    if length(response_times)>0
        T_set_time = Image_Time(response_times(1));
    elseif length(response_times)==0
        T_set_time = 1;
    end
    
    %I am hoping I can ignore all this YFP stuff and go straight to the
    %else
    if exist('clean_YFP')==1
        CFPdivYFP = clean_CFP(1:20)./clean_YFP(1:20);%%%this is just for visualization....
        mean_factor = mean(CFPdivYFP);                 %%%Corrects the smaller signal to
        YFP_corrected = clean_YFP*mean_factor;        %%%be of equal beginnning amp as the larger one.
        diff_YFP = diff(YFP_corrected);% All of the YFP_corrected stuff won't generate numbers
        diff_CFP = diff(clean_CFP);   % in the end, its just for display.
        slope_diff = [];
        for i = 1:length(diff_CFP)
            slope_diff(end+1) = diff_CFP(i)-diff_YFP(i);
        end
        
        figure();
        set(0,'DefaultFigurePosition',[100 100 2000 2000]);
        subplot(2,2,1)
        plotyy(Image_Time, clean_CFP, Image_Time, YFP_corrected)
        yL = get(gca,'YLim');
        line([T_set_time T_set_time],yL,'Color','r'); %Red line is a reference for the first fret peak
        title('CFP and YFP traces')
        
        subplot(2,2,2)
        plotyy(Image_Time(1:end-1), diff_CFP, Image_Time(1:end-1), diff_YFP)
        yyL = get(gca,'YLim');
        line([T_set_time T_set_time],yyL,'Color','r');
        title('Slopes of CFP/YFP. Find the first reciprocal change, and double click. Red line indicates position of first fret peak.')
        hr = gca;
        
        subplot(2,2,3)
        plot(slope_diff)
        yyyL = get(gca,'YLim');
        line([T_set_time T_set_time],yyyL,'Color','r');
        title('Difference in slope between CFP and YFP')
        
        subplot(2,2,4)
        plotyy(Image_Time, del_F, Image_Time, raw_temperature)
        hold on
        plot(dff,'r')
        yyyyL = get(gca,'YLim');
        line([T_set_time T_set_time],yyyyL,'Color','r');
        title('Fret Signal')
        [recip_x, recip_y] = getpts(hr);
    else
        figure('position', [1200 0 700 1000]);
        plotyy(Image_Time, del_F, Image_Time, raw_temperature)
        hold on
        plot(dff,'r');%ylim([-20 inf]);
        yyyyL = get(gca,'YLim');
        line([T_set_time T_set_time],yyyyL,'Color','r');
        title('Double click where the response begins')
        hr = gca;
        [recip_x, recip_y] = getpts(hr);
    end
    %I guess the output of this ends up being recip_x&y
    
    
    recip_x = round(recip_x);
    close
    
    if exist('BLC_clean_delta_F')==1
        plot(BLC_raw_delta_F);
        title('Click on the min, then the max of the signal peak, then double click white space.') 
        h = gca;
    else
        plot(dff);
        title('Click on the min, then the max of the signal peak, then double click white space.') 
        h = gca;
    end
    %Not sure why this is needed, output is h which is the points you
    %clicked on
    
    %should be able to get rid of all the c iteration stuff
    [locs, pks] = getpts(h);
    fret.pks = pks(2)-pks(1);
    fret.locs = locs(1);
    close
    
    %can also get rid of end+1 stuff
    T_star_time = round((recip_x));
    T_star_temp = raw_temperature(T_star_time);
                                                      
    T_set_time = round(T_set_time);                                                
    
    
%     %is this necessary??
%     if T_set_time>0
%         T_set_point(end+1) = raw_temperature(T_set_time);% correcting this once we have
%     elseif T_set_time<=0                                 % mean_peak_diff
%         T_set_point(end+1) =nan;
%     end
    
    %commenting this out cuz not iterating
%     T_star_temps(c)
%     T_set_point(c)
%     c = c + 1;

total_fret_amp = [fret.pks];

end