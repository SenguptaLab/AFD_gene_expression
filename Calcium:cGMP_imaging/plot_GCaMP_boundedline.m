%for plotting averages of all dff traces in a condition as well as
%standard error shading
%requires https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m

%I keep my aggregated deltaF/F traces (the variable dff from
%GCaMP_vs_temperature_NH) in excel spreadsheets. To bring them here I just
%make variables for each condition - dff_A = []; dff_B = []; and so on and
%paste in the traces

%average the dffs
average_trace_A = mean(dff_A,2);
average_trace_B = mean(dff_B,2);
average_trace_C = mean(dff_C,2);
average_trace_D = mean(dff_D,2);    

%get standard deviation at every frame
for i=1:length(dff_A)
    std_A(i) = std(dff_A(i,:));
end

for i=1:length(dff_B)
    std_B(i) = std(dff_B(i,:));
end

for i=1:length(dff_C)
    std_C(i) = std(dff_C(i,:));
end

for i=1:length(dff_D)
     std_D(i) = std(dff_D(i,:));
end

%calculate standard errors
sem_A = std_A/sqrt(size(dff_A,2));
sem_B = std_B/sqrt(size(dff_B,2));
sem_C = std_C/sqrt(size(dff_C,2));
sem_D = std_D/sqrt(size(dff_D,2));

%plot average traces, SEMs, and temperature ramp (use real_temp_ramp
%variable from GCaMP_vs_temperature_NH)
%adjust ylim as desired
%also adjust length of plotting according to your video length (typically
%1:141 for mine)
figure('position', [1200 0 700 1000]); hold on;
yyaxis left; ylim([-5, 30]);
ylabel('\deltaF/F (%)');
boundedline((1:141),average_trace_A(1:141),sem_A(1:141),'k','alpha','transparency',.1);
boundedline((1:141),average_trace_B(1:141),sem_B(1:141),'r','alpha','transparency',.1);
boundedline((1:141),average_trace_C(1:141),sem_C(1:141),'b','alpha','transparency',.1);
boundedline((1:141),average_trace_D(1:141),sem_D(1:141),'g','alpha','transparency',.1);

yyaxis right; ylim([(min(real_temp)-1), (max(real_temp)+1)]); 
plot(real_temp_ramp(1:141),'-g','LineWidth',2);
ylabel('Temperature');

