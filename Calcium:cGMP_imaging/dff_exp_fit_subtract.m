%fits an exponential to the regions of dff with no signal (that
%just have photobleaching), then subtracts that fit. 
%Probably only necessary for FlincG3 imaging, not GCaMP
%Used for beautifying traces before making figure, not for calculating
%T* response threshold
dff_no_start = dff(21:length(dff),:);
dff_sub = [];

%get frame#s before and after signal (to exclude the response from the
%curve fit
pre_response = input('time before response starts?');
if pre_response < 21
    pre_response = 21;
end
post_response = input('time after response starts?');
for i = 1:size(dff_no_start,2)
    skip_times = [(pre_response-20):(post_response-20)];
    dff_skip = dff_no_start(:,i);
    dff_skip(skip_times)=NaN;

    %matching times
    t = [1:(length(dff)-20)];

    %get indices for parts of curve with data points, these will be used
    %for the fit
    idxValid = ~isnan(dff_skip);

    %do the fit and extract the curve over t
    dff_fit_function = fit(t(idxValid).',dff_skip(idxValid),'exp2');
    dff_fit_curve = feval(dff_fit_function,t.');
    
    %subtract curve fit
    dff_sub = [dff_sub, (dff_no_start(:,i)-dff_fit_curve)];
end

average_trace = mean(dff_sub,2);

temp_ramp_no_start = real_temp_ramp(21:length(real_temp_ramp));

figure('position', [1200 0 700 1000]); hold on;

%subplot(size(rois,2),1,i);
yyaxis left; ylim([-10, 20]);%ylim([min(dff_sub(:)), (max(dff_sub(:)))]); 
F_trace = plot(dff_sub,'-k'); %F_trace.Color(4)=0.3;
%F2_trace = plot(dff(:,2),'-r'); F2_trace.Color(4)=0.3;
plot(average_trace,'-k','LineWidth',2);
ylabel('\deltaF/F (%)');
yyaxis right; ylim([(min(temp_ramp_no_start)-1), (max(temp_ramp_no_start)+1)]); 
plot(temp_ramp_no_start,'-g','LineWidth',2);
%or use this for average-
%yyaxis right; ylim([(min(average_T)-1), (max(average_T)+1)]); 
%plot(average_T,'-g','LineWidth',2);
ylabel('Temperature');


delta_F_exp_corrected = strcat(folder, 'delta_F_exp_corrected');
saveas(gcf, delta_F_exp_corrected,'fig');