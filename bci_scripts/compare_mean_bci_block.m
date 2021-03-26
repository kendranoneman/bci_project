clear;
clear figs;

%monkey = 'Pepe';
monkey = 'Satchel';

load(sprintf('%s/bcidistance.mat',monkey));
load(sprintf('%s/bciacqtime.mat',monkey)); 
load(sprintf('%s/blockshamacqtime.mat',monkey));
load(sprintf('%s/blockshamdistance.mat',monkey));
%load(sprintf('%s/bcishamacqtime.mat',monkey));
%load(sprintf('%s/bcishamdistance.mat',monkey));
load(sprintf('%s/sessindindex.mat',monkey));

unique_sessidx = unique(sessindindex);
if not(isfolder('hist_figs'))
    mkdir('hist_figs')
end

unique_sessidx = unique(sessindindex);

%% Compute the difference in mean between BCI and block sham
% Store these in a vector that is 1 by the number of sessions
correct_diff_all = [];
error_diff_all = [];
tot_diff_all = [];

bci_tot_acq = [];
block_tot_acq = [];
bci_tot_dists_all = [];
block_tot_dists_all = [];
for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);
    
    bci_acqtimes = bciacqtime(sessindindex==s,:); 
    bci_acqtimes = bci_acqtimes(:,1:80);
    block_acqtimes = blockshamacqtime(sessindindex==s,:); 
    
    bci_tot_acq = [bci_tot_acq; bci_acqtimes];
    block_tot_acq = [block_tot_acq; block_acqtimes];
    
    bci_correct_idx = ~isnan(bci_acqtimes);
    block_correct_idx = ~isnan(block_acqtimes);
    
    bcidists = bcidistance{i};
    bcidists = bcidists(:,1:80);
    blockdists = blockshamdistance{i}; 
    
    bci_trim = [];
    for q2 = 1:size(bcidists,1)
        for m2 = 1:size(bcidists,2)
            aa2 = bci_acqtimes(q2,m2);
            bb2 = bcidists(q2,m2);
            if isnan(aa2)
                bci_trim_dist2=nan;
            else
                mat2 = bb2{1}([1:aa2]);
                bci_trim_dist2 = mat2cell(mat2,1,aa2);
                bci_trim = [bci_trim; bci_trim_dist2];
            end
        end
    end
   
    block_trim = [];
    for q = 1:size(blockdists,1)
        for m = 1:size(blockdists,2)
            aa = block_acqtimes(q,m);
            bb = blockdists(q,m);
            if isnan(aa)
                block_trim_dist=nan;
            else
                mat = bb{1}([1:aa]);
                block_trim_dist = mat2cell(mat,1,aa);
                block_trim = [block_trim; block_trim_dist];
            end
        end
    end

    % all trials
    bci_correct_dists = bci_trim;
    bci_error_dists = bcidists(~bci_correct_idx);
    bci_tot_dists = [bci_correct_dists; bci_error_dists];
    block_correct_dists = block_trim;
    block_error_dists = blockdists(~block_correct_idx);
    block_tot_dists = [block_correct_dists; block_error_dists];

    % mean of each trial
    bci_correct_all_mean = [];
    for j = 1:length(bci_correct_dists)
        bci_correct_all_mean = [bci_correct_all_mean, mean(bci_correct_dists{j})];
    end
    bci_correct_mean = mean(bci_correct_all_mean);

    bci_error_all_mean = [];
    for a = 1:length(bci_error_dists)
        bci_error_all_mean = [bci_error_all_mean, mean(bci_error_dists{a})];
    end
    bci_error_mean = mean(bci_error_all_mean);
    
    bci_tot_all_mean = [];
    for q = 1:length(bci_tot_dists)
        bci_tot_all_mean = [bci_tot_all_mean, mean(bci_tot_dists{q})];
    end
    bci_tot_dists_all = [bci_tot_dists_all , bci_tot_all_mean];

    block_correct_all_mean = [];
    for l = 1:length(block_correct_dists)
        block_correct_all_mean = [block_correct_all_mean, mean(block_correct_dists{l})];
    end
    block_correct_mean = mean(block_correct_all_mean);

    block_error_all_mean = [];
    for c = 1:length(block_error_dists)
        block_error_all_mean = [block_error_all_mean, mean(block_error_dists{c})];
    end
    block_error_mean = mean(block_error_all_mean);
    
    block_tot_all_mean = [];
    for p = 1:length(block_tot_dists)
        block_tot_all_mean = [block_tot_all_mean, mean(block_tot_dists{p})];
    end
    block_tot_dists_all = [block_tot_dists_all , block_tot_all_mean];

    correct_diff = block_correct_mean - bci_correct_mean;
    error_diff = block_error_mean - bci_error_mean;
    %tot_diff = block_tot_mean - bci_tot_mean;

    correct_diff_all = [correct_diff_all ; correct_diff];
    error_diff_all = [error_diff_all ; error_diff];
    %tot_diff_all = [tot_diff_all ; tot_diff];
    
%     bci_tot_dists_all = [bci_tot_dists_all ; bci_tot_dists];
%     block_tot_dists_all = [block_tot_dists_all ; block_tot_dists];
end

correct_diff_mean = mean(correct_diff_all,'omitnan');
error_diff_mean = mean(error_diff_all,'omitnan');
%tot_diff_mean = mean(tot_diff_all, 'omitnan');

bci_correct_acq_mean = mean(bci_tot_acq, 'omitnan');
block_correct_acq_mean = mean(block_tot_acq, 'omitnan');
%% Plot that vector as a histogram 
% Do statistical test to see if histogram is significantly different from 0

fontsize = 16;

binwidth = 0.02;
blue = [0, 0.4470, 0.7410];
orange = [0.8500 0.3250 0.0980];
green = [0.4660, 0.6740, 0.1880];	
gray = [0.17 0.17 0.17];

h=figure('Visible','Off');
hold on
histogram(correct_diff_all,'BinWidth',binwidth,'FaceAlpha',0.5,'EdgeColor',blue, 'LineWidth', 1.25);
histogram(error_diff_all,'BinWidth',binwidth,'FaceAlpha',0.5,'EdgeColor',orange,'LineWidth', 1.25);
xline(correct_diff_mean,'Color',blue,'LineWidth',2);
xline(error_diff_mean,'Color',orange,'LineWidth',2);
xline(0,'k--','LineWidth',1.5);
%title(sprintf('Difference in Mean Distance between BCI & Block Sham Across All Sessions: %s',monkey));
xlabel('Difference in Mean Distance (Block Sham - BCI)','Fontsize',fontsize)
ylabel('# of Sessions','Fontsize',fontsize)
xlim([-0.1 0.13])
ylim([0 15])
legend('Correct','Missed');
saveas(h,sprintf('hist_figs/%s/hist-bci-block-diff.jpg',monkey));
hold off

% h2=figure('Visible','Off');
% hold on
% histogram(tot_diff_all, 'BinWidth',binwidth,'FaceColor',gray,'FaceAlpha',0.1,'EdgeColor','k','LineWidth', 1.25)
% xline(tot_diff_mean, 'Color','k', 'LineWidth',2);
% xline(0,'k--','LineWidth',1.5);
% %title(sprintf('Difference in Mean Distance between BCI & Block Sham Across All Sessions: %s',monkey));
% xlabel('Difference in Mean Distance (Block Sham - BCI)','Fontsize',fontsize)
% ylabel('# of Sessions','Fontsize',fontsize)
% xlim([-0.1 0.13])
% ylim([0 15])
% legend('All');
% saveas(h2,sprintf('hist_figs/%s/hist-bci-block-diff-tot.jpg',monkey));

%% Histogram of Acquisition Times
% h3 = figure('Visible','Off');
% hold on
% binwidth = 0.02;
% histogram(bci_correct_acqtimes,'k','BinWidth',binwidth, 'FaceAlpha',0.5,'LineWidth', 1.25);
% histogram(block_correct_acqtimes, 'm','BinWidth',binwidth,'FaceAlpha',0.5,'LineWidth', 1.25);
% xline(bci_correct_acq_mean, 'Color','k', 'LineWidth',2);
% xline(block_correct_acq_mean, 'Color','m', 'LineWidth',2);
% xline(0,'k--','LineWidth',1.5);
% title(sprintf('Acquisition Times : %s',monkey));
% xlabel('Difference in Mean Distance (Block Sham - BCI)')
% ylabel('# of Sessions')
% xlim([-0.03 0.13])
% ylim([0 12])
% legend('All');
% saveas(h2,sprintf('hist_figs/%s/hist-bci-block-diff-tot.jpg',monkey));
