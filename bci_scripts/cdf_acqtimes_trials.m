clear;
clear figs;

load('bciacqtime.mat'); %trials vs. blocks
load('bcishamacqtime'); 
load('blockshamacqtime');
load('sessindindex.mat');

unique_sessidx = unique(sessindindex);
if not(isfolder('cdf_figs'))
    mkdir('cdf_figs')
end

for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);

    %Acquisition times for each session
    bci_acqtimes = bciacqtime(sessindindex==s,:);
    sham_acqtimes = bcishamacqtime(sessindindex==s,:);
    block_acqtimes = blockshamacqtime(sessindindex==s,:);

    %Tossing out the missed trials (bin = 50 ms)
    bci_correct = (bci_acqtimes(~isnan(bci_acqtimes))) * 50; %ms
    sham_correct = (sham_acqtimes(~isnan(sham_acqtimes))) * 50;
    block_correct = (block_acqtimes(~isnan(block_acqtimes))) * 50;cd

    %Calculating ecdf
    [f1_bci,x1] = ecdf(bci_correct);
    [f1_sham,x2] = ecdf(sham_correct);
    [f1_block,x3] = ecdf(block_correct);

    %Plotting ecdf
    h=figure('Visible','Off');
    hold on 
    plot(x1,f1_bci*90,'r');
    plot(x2,f1_sham*10, 'b');
    plot(x3,f1_block*20, 'm');
    title(sprintf('ECDF of Acquisition Times (session = %d)', s));
    xlabel('x (ms)');
    ylabel('F(x) * (# of Trials)');
    legend('BCI (90 trials)','BCI Sham (10 trials)','Block Sham (20 trials)','Location','Best');
    saveas(h,sprintf('cdf_figs/ecdf-trials-s%d.png',s));
end