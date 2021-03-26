clear;
clear figs;

%monkey = 'Pepe';
monkey = 'Satchel';

load(sprintf('%s/bcidistance.mat',monkey));
load(sprintf('%s/bciacqtime.mat',monkey)); 
load(sprintf('%s/blockshamacqtime.mat',monkey));
load(sprintf('%s/blockshamdistance.mat',monkey));
load(sprintf('%s/sessindindex.mat',monkey));

unique_sessidx = unique(sessindindex);
if not(isfolder('hist_figs'))
    mkdir('hist_figs')
end

unique_sessidx = unique(sessindindex);

%% Acquisition Times
bci_tot_acq = [];
block_tot_acq = [];
bci_sem_tot = [];
block_sem_tot = [];

for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);
    
    bci_acqtimes = bciacqtime(sessindindex==s,:); 
    bci_acqtimes = bci_acqtimes(:,1:80);
    block_acqtimes = blockshamacqtime(sessindindex==s,:); 
    
    bci_tot_acq = [bci_tot_acq; bci_acqtimes];
    block_tot_acq = [block_tot_acq; block_acqtimes];
end

bci_correct_acq_mean = mean(bci_tot_acq, 2, 'omitnan');
block_correct_acq_mean = mean(block_tot_acq, 2, 'omitnan');

bci_mean = mean(bci_correct_acq_mean,'omitnan');
block_mean = mean(block_correct_acq_mean,'omitnan');

%% Histogram of Acquisition Times
h3 = figure;
hold on
binwidth = 2;
histogram(bci_correct_acq_mean,'FaceColor','k','FaceAlpha',0.2,'BinWidth',binwidth,'LineWidth', 1.25);
histogram(block_correct_acq_mean,'FaceColor', 'm','FaceAlpha',0.2,'BinWidth',binwidth,'LineWidth', 1.25);
xline(bci_mean, 'k', 'LineWidth',2);
xline(block_mean, 'm', 'LineWidth',2);
%xline(0,'k--','LineWidth',1.5);
title(sprintf('Acquisition Times for Correct Trials Across All Sessions: %s',monkey));
xlabel('Acquisition Times')
ylabel('# of Trials')
xlim([0 70])
ylim([0 55])
legend('BCI','Block Sham');
saveas(h3,sprintf('hist_figs/%s/hist-acq.jpg',monkey));

%% Acquisition Times over course of Trial
% x = (1:75)*50;
% plot(x,bci_correct_acq_mean,'k');
% plot(x,block_correct_acq_mean, 'm');
