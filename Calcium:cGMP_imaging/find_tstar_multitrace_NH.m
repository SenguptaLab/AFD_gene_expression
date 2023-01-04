%This function feeds all your deltaF/F traces in dff variable through
%T_star_loop, outputting a list of TStars
%for each trace, double click where T star is, then click the bottom
%and top of the trace to get amplitude (which we don't use anyway)

T_star_temps = [];
total_dff_amps = [];

for i = 1:length(dff(1,:))
    
    T_star_temp = [];
    
    Image_Time = 1:length(real_temp_ramp);
    raw_temperature = real_temp_ramp;

    steep_points = [];             % This is the important part for finding the first 
    time_steep = [];               % sustained increase in the del_F signal.
    increasing = [];
    del_F = [];
    diff_del_F = [];
    
    %see function T_star_loop
    [T_star_temp,total_dff_amp] = T_star_loop(Image_Time,raw_temperature,dff(:,i));
    
    T_star_temps(i) = T_star_temp;
    total_dff_amps(i) = total_dff_amp;
    
    
end
 
T_star_WT = T_star_temps(1);
T_star_mut = T_star_temps(2);

average_t_star = mean(T_star_temps)
save analysis