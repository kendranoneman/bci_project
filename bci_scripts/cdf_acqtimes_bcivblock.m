clear;
clear figs;

%% Specifying monkey, loading data, formatting into correct/missed
monkey = 'Pepe';
%monkey = 'Satchel';

load(sprintf('%s/bcidistance.mat',monkey));
load(sprintf('%s/bciacqtime.mat',monkey)); 
load(sprintf('%s/blockshamacqtime.mat',monkey));
load(sprintf('%s/blockshamdistance.mat',monkey));
load(sprintf('%s/bcishamacqtime.mat',monkey));
load(sprintf('%s/bcishamdistance.mat',monkey));
load(sprintf('%s/sessindindex.mat',monkey));

% Makes a new folder for these if necessary
unique_sessidx = unique(sessindindex);
if not(isfolder('cdf_figs'))
    mkdir('cdf_figs')
end

bci_all = [];
sham_all = [];
block_all = [];
for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);

    % Acquisition times for each session
    bci_acqtimes = bciacqtime(sessindindex==s,:);
    sham_acqtimes = bcishamacqtime(sessindindex==s,:);
    block_acqtimes = blockshamacqtime(sessindindex==s,:);

    % Tossing out the missed trials (bin = 50 ms)
    bci_correct = (bci_acqtimes(~isnan(bci_acqtimes))) * 50; %ms
    sham_correct = (sham_acqtimes(~isnan(sham_acqtimes))) * 50;
    block_correct = (block_acqtimes(~isnan(block_acqtimes))) * 50;
    
    bci_all = [bci_all; bci_correct]; 
    sham_all = [sham_all; sham_correct]; 
    block_all = [block_all; block_correct]; 
end

bci_f_all = []; block_f_all = [];
bci_ci_all = []; block_ci_all = [];

%% Randomly sampling from BCI data
% BCI = 90 blocks, Sham = 10 blocks, Block = 20 blocks

n=20; % number of times to sample bci/block
for i = 1:n
    bci_samp = datasample(bci_all,length(sham_all)); 
    block_samp = datasample(block_all, length(sham_all)); 
    
    [f_bci,x1] = ecdf(bci_samp);
    bci_f_all = [bci_f_all , f_bci]; 
    [f_block,x2] = ecdf(block_samp); 
    block_f_all = [block_f_all, f_block]; 
end

[f_sham,x3] = ecdf(sham_all);
bci_f_mean = mean(bci_f_all,2);
block_f_mean = mean(block_f_all,2);

% Calculating 5th/95th percentiles
bci_f_ci = prctile(bci_f_all,[5 95],2);
block_f_ci = prctile(block_f_all,[5 95],2);


%% Plotting empirical CDFs for all three trial types
% Plotting ecdf

fontsize = 16;

h=figure('Visible','Off');
hold on 
plot(x1,bci_f_mean, 'k','LineWidth',1.5);
plot(x3,f_sham, 'b', 'LineWidth',1.5);
plot(x2,block_f_mean, 'm','LineWidth',1.5);
plot(x1,bci_f_ci(:,1), 'k--');
plot(x1,bci_f_ci(:,2), 'k--');
plot(x2,block_f_ci(:,1), 'm--');
plot(x2,block_f_ci(:,2), 'm--');
%title(sprintf('ECDF of Acquisition Times w/ 5th & 95th Percentiles: %s',monkey));
xlabel('x (ms)','FontSize',fontsize);
ylabel('F(x)','FontSize',fontsize);
xlim([500 4000]);
ylim([0 1]);
legend('BCI','BCI Sham','Block Sham','Location','southeast');
saveas(h,sprintf('cdf_figs/%s/ecdf-confint.png',monkey));

% Plotting zoomed in view of ecdf 
j=figure('Visible','Off');
hold on 
plot(x1(20:40),bci_f_mean(20:40), 'k','LineWidth',1.5);
plot(x3(20:40),f_sham(20:40), 'b', 'LineWidth',1.5);
plot(x2(20:40),block_f_mean(20:40), 'm','LineWidth',1.5);
plot(x1(20:40),bci_f_ci(20:40,1), 'k--');
plot(x1(20:40),bci_f_ci(20:40,2), 'k--');
plot(x2(20:40),block_f_ci(20:40,1), 'm--');
plot(x2(20:40),block_f_ci(20:40,2), 'm--');
%title(sprintf('ECDF of Acquisition Times w/ 5th & 95th Percentiles: %s',monkey));
xlabel('x (ms)','FontSize',fontsize);
ylabel('F(x)','FontSize',fontsize);
xlim([1600 2800]);
ylim([0.4 0.85]);
legend('BCI','BCI Sham','Block Sham','Location','southeast');
saveas(j,sprintf('cdf_figs/%s/ecdf-confint-zoom.png',monkey));