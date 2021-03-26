clear;
clear figs;

load('bcidistance.mat');
load('bciacqtime.mat'); 
load('blockshamacqtime.mat');
load('blockshamdistance.mat');
load('bcishamacqtime.mat')
load('bcishamdistance.mat')
load('sessindindex.mat');

unique_sessidx = unique(sessindindex);
if not(isfolder('hist_figs'))
    mkdir('hist_figs')
end

unique_sessidx = unique(sessindindex);

s=3;
bci_acqtimes = bciacqtime(sessindindex==s,:); 
bcidists = bcidistance{i};
block_acqtimes = blockshamacqtime(sessindindex==s,:); 
blockdists = blockshamdistance{i}; 
sham_acqtimes = bcishamacqtime(sessindindex==s,:); 
shamdists = bcishamdistance{i}; 

bci_correct_idx = ~isnan(bci_acqtimes);
block_correct_idx = ~isnan(block_acqtimes);
sham_correct_idx = ~isnan(sham_acqtimes);

bci_correct_dists = bcidists(bci_correct_idx);
bci_error_dists = bcidists(~bci_correct_idx);
sham_correct_dists = bcidists(sham_correct_idx);
sham_error_dists = bcidists(~sham_correct_idx);
block_correct_dists = bcidists(block_correct_idx);
block_error_dists = bcidists(~block_correct_idx);

% iterate and concatenate all the vectors for correct vs error into a
% single histogram and plot the others
bci_correct_all = [];
for i = 1:length(bci_correct_dists)
    bci_correct_all = [bci_correct_all , bci_correct_dists{i}];
end
bci_error_all = [];
for j = 1:length(bci_error_dists)
    bci_error_all = [bci_error_all , bci_error_dists{j}];
end

sham_correct_all = [];
for k = 1:length(sham_correct_dists)
    sham_correct_all = [sham_correct_all , sham_correct_dists{k}];
end
sham_error_all = [];
for l = 1:length(sham_error_dists)
    sham_error_all = [sham_error_all , sham_error_dists{l}];
end

block_correct_all = [];
for p = 1:length(block_correct_dists)
    block_correct_all = [block_correct_all , block_correct_dists{p}];
end
block_error_all = [];
for q = 1:length(sham_error_dists)
    block_error_all = [block_error_all , block_error_dists{q}];
end

h=figure('Visible','off');
hold on
histogram(bci_correct_all);
histogram(bci_error_all);
title(sprintf('BCI: Correct vs. Missed Distances Over All Trials (session = %d)', s));
legend('Correct','Missed','Location','Best');
saveas(h,'hist_figs/hist-bci-s3.jpg');

p=figure('Visible','off');
hold on
histogram(sham_correct_all);
histogram(sham_error_all);
title(sprintf('BCI Sham: Correct vs. Missed Distances Over All Trials (session = %d)', s));
legend('Correct','Missed','Location','Best');
saveas(p,'hist_figs/hist-sham-s3.jpg');

q=figure('Visible','off');
hold on
histogram(block_correct_all);
histogram(block_error_all);
title(sprintf('Block Sham: Correct vs. Missed Distances Over All Trials (session = %d)', s));
legend('Correct','Missed','Location','Best');
saveas(q,'hist_figs/hist-block-s3.jpg');

%% compare distances on all trials across bci vs bci sham vs block sham (bci < bci sham < block sham)
bci_all = [bci_correct_all , bci_error_all]; %30561
sham_all = [sham_correct_all, sham_error_all]; %3910
block_all = [block_correct_all, block_error_all]; %4586

a = figure;
hold on
histogram(bci_all);
histogram(sham_all);
histogram(block_all);
title(sprintf('Distances on all trials across BCI vs. BCI Sham vs. Block Sham (session = %d)',s));
legend('BCI','BCI Sham', 'Block Sham');
saveas(a,'hist_figs/hist-all-s3.jpg');

%% compare distances on correct across bci vs bci sham vs block sham
b = figure;
hold on
histogram(bci_correct_all);
histogram(sham_correct_all);
histogram(block_correct_all);
title(sprintf('Distances all correct trials across BCI vs. BCI Sham vs. Block Sham (session = %d)',s));
legend('BCI', 'BCI Sham', 'Block Sham');
saveas(b,'hist_figs/hist-correct-s3.jpg');

%% compare distances on error across bci vs bci sham vs block sham
c = figure;
hold on
histogram(bci_error_all);
histogram(sham_error_all);
histogram(block_error_all);
title(sprintf('Distances all missed trials across BCI vs. BCI Sham vs. Block Sham (session = %d)',s));
legend('BCI', 'BCI Sham', 'Block Sham');
saveas(b,'hist_figs/hist-error-s3.jpg');
