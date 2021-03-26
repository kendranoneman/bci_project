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

v = VideoWriter ('cdf_figs/test.avi');
v.FrameRate = 3;
open(v);
h=figure(); hold on
for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);

    %Acquisition times for each session
    bci_acqtimes = bciacqtime(sessindindex==s,:);
    sham_acqtimes = bcishamacqtime(sessindindex==s,:);
    block_acqtimes = blockshamacqtime(sessindindex==s,:);

    %Tossing out the missed trials (bin = 50 ms)
    bci_correct = (bci_acqtimes(~isnan(bci_acqtimes))) * 50; %ms
    sham_correct = (sham_acqtimes(~isnan(sham_acqtimes))) * 50;
    block_correct = (block_acqtimes(~isnan(block_acqtimes))) * 50;

    %Calculating ecdf
    [f_bci,x1] = ecdf(bci_correct);
    [f_sham,x2] = ecdf(sham_correct);
    [f_block,x3] = ecdf(block_correct);

    %Plotting ecdf
    plot(x1,f_bci,'k');
    plot(x2,f_sham, 'b');
    plot(x3,f_block, 'g');
    title('ECDF of Acquisition Times');
    xlabel('x (ms)');
    ylabel('F(x)');
    legend('BCI','BCI Sham','Block Sham','Location','Best');
    %saveas(h,'cdf_figs/ecdf-stack.png');
    
    frame = getframe (gcf);
    writeVideo (v, frame);
end
close(v);