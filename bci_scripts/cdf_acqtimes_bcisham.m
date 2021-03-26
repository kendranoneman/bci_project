clear;
clear figs;

load('bcishamacqtime'); 
load('sessindindex.mat');

unique_sessidx = unique(sessindindex);
if not(isfolder('cdf_figs'))
    mkdir('cdf_figs')
end

sham_all = [];
for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);

    %Acquisition times for each session
    sham_acqtimes = bcishamacqtime(sessindindex==s,:);

    %Tossing out the missed trials (bin = 50 ms)
    sham_correct = (sham_acqtimes(~isnan(sham_acqtimes))) * 50;
    sham_all = [sham_all; sham_correct]
    
    %Calculating ecdf
    [f_sham,x] = ecdf(sham_all);
end

%Plotting ecdf
h=figure('Visible','Off');
hold on 
plot(x,f_sham, 'b');
title('ECDF of Acquisition Times: BCI Sham');
xlabel('x (ms)');
ylabel('F(x)');
legend('BCI Sham','Location','Best');
saveas(h,'cdf_figs/ecdf-bcisham.png');