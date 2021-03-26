clear;
clear figs;

%% Selecting monkey, loading data, separating into correct/missed
monkey = 'Pepe';
%monkey = 'Satchel';

load(sprintf('%s/bcidistance.mat',monkey));
load(sprintf('%s/bciacqtime.mat',monkey)); 
load(sprintf('%s/blockshamacqtime.mat',monkey));
load(sprintf('%s/blockshamdistance.mat',monkey));
load(sprintf('%s/bcishamacqtime.mat',monkey));
load(sprintf('%s/bcishamdistance.mat',monkey));
load(sprintf('%s/sessindindex.mat',monkey));

% Create folder if doesn't exist
unique_sessidx = unique(sessindindex);
if not(isfolder('hist_figs'))
    mkdir('hist_figs')
end

% Using sessindindex file to determine session numbers
unique_sessidx = unique(sessindindex);

bci_correct_all = [];
bci_error_all = [];
sham_correct_all = [];
sham_error_all = [];
block_correct_all = [];
block_error_all = [];


for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);
    
    bcidists = bcidistance{i};
    bcidists = bcidists(:,1:80);
    blockdists = blockshamdistance{i};
    shamdists = bcishamdistance{i}; 
    
    bci_acqtimes = bciacqtime(sessindindex==s,:); 
    bci_acqtimes = bci_acqtimes(:,1:80); 
    block_acqtimes = blockshamacqtime(sessindindex==s,:); 
    sham_acqtimes = bcishamacqtime(sessindindex==s,:);

    bci_correct_idx = ~isnan(bci_acqtimes);
    sham_correct_idx = ~isnan(sham_acqtimes);
    block_correct_idx = ~isnan(block_acqtimes);
    
    % Matching up size of acquisition times & distances
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
    
    sham_trim = [];
    for q3 = 1:size(shamdists,1)
        for m3 = 1:size(shamdists,2)
            aa3 = sham_acqtimes(q3,m3);
            bb3 = shamdists(q3,m3);
            if isnan(aa3)
                sham_trim_dist2=nan;
            else
                mat2 = bb3{1}([1:aa3]);
                sham_trim_dist = mat2cell(mat2,1,aa3);
                sham_trim = [sham_trim; sham_trim_dist];
            end
        end
    end
    
    % Most important for block trim, since the trial is made up of 
    % completely sham trals; doesn't stop when monkey gets correct
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
    sham_correct_dists = sham_trim;
    sham_error_dists = shamdists(~sham_correct_idx);
    block_correct_dists = block_trim;
    block_error_dists = blockdists(~block_correct_idx);
     

    % mean of each trial
    bci_correct_all_mean = [];
    for j = 1:length(bci_correct_dists)
        bci_correct_all_mean = [bci_correct_all_mean, mean(bci_correct_dists{j})];
    end
    bci_correct_mean = mean(bci_correct_all_mean,'omitnan');
    bci_correct_all = [bci_correct_all, bci_correct_all_mean];

    bci_error_all_mean = [];
    for a = 1:length(bci_error_dists)
        bci_error_all_mean = [bci_error_all_mean, mean(bci_error_dists{a})];
    end
    bci_error_mean = mean(bci_error_all_mean,'omitnan');
    bci_error_all = [bci_error_all , bci_error_all_mean];

    sham_correct_all_mean = [];
    for k = 1:length(sham_correct_dists)
        sham_correct_all_mean = [sham_correct_all_mean, mean(sham_correct_dists{k})];
    end
    sham_correct_mean = mean(sham_correct_all_mean,'omitnan');
    sham_correct_all = [sham_correct_all, sham_correct_all_mean];

    sham_error_all_mean = [];
    for b = 1:length(sham_error_dists)
        sham_error_all_mean = [sham_error_all_mean, mean(sham_error_dists{b})];
    end
    sham_error_mean = mean(sham_error_all_mean,'omitnan');
    sham_error_all = [sham_error_all, sham_error_all_mean];

    block_correct_all_mean = [];
    for l = 1:length(block_correct_dists)
        block_correct_all_mean = [block_correct_all_mean, mean(block_correct_dists{l})];
    end
    block_correct_mean = mean(block_correct_all_mean,'omitnan');
    block_correct_all = [block_correct_all, block_correct_all_mean];

    block_error_all_mean = [];
    for c = 1:length(block_error_dists)
        block_error_all_mean = [block_error_all_mean, mean(block_error_dists{c})];
    end
    block_error_mean = mean(block_error_all_mean,'omitnan');
    block_error_all = [block_error_all, block_error_all_mean];
    
    %% compare distances on error across bci vs block sham
    % for each session
    
    fontsize = 16;
    
    if not(isfolder(sprintf('hist_figs/%s/bci_block-error',monkey)))
        mkdir(sprintf('hist_figs/%s/bci_block-error',monkey))
    end
    
    binwidth = 0.08;
    blue = [0, 0.4470, 0.7410];
    orange = [0.8500 0.3250 0.0980];
    green = [0.4660, 0.6740, 0.1880];	
    gray = [0.17 0.17 0.17];
    
    f = figure('Visible','off');
    hold on
    h1 = histogram(bci_error_all_mean,'BinWidth',binwidth,'FaceColor','k');
    h2 = histogram(block_error_all_mean,'BinWidth',binwidth,'FaceColor','m');
    xline(bci_error_mean,'k','LineWidth',2);
    xline(block_error_mean,'m','LineWidth',2);
    %title(sprintf('Average Distances all missed trials across BCI vs. Block Sham: %s (session = %d)',monkey, s));
    xlabel('Distance from Target','FontSize',fontsize)
    ylabel('# of Trials','FontSize',fontsize)
    xlim([0.8 1.93])
    ylim([0 210])
    legend('BCI','Block Sham',sprintf('Mean = %.3f',bci_error_mean), sprintf('Mean = %.3f',block_error_mean));
    saveas(f,sprintf('hist_figs/%s/bci_block-error/hist-error-mean-s%d.jpg',monkey,s));
    
    %% compare distances on correct across bci vs block sham
    % individual sessions
    if not(isfolder(sprintf('hist_figs/%s/bci_block-correct',monkey)))
        mkdir(sprintf('hist_figs/%s/bci_block-correct',monkey))
    end
    
    e = figure('Visible','off');
    hold on
    h3 = histogram(bci_correct_all_mean,'BinWidth',binwidth,'FaceColor','k');
    h4 = histogram(block_correct_all_mean,'BinWidth',binwidth,'FaceColor','m');
    xline(bci_correct_mean,'k','LineWidth',2);
    xline(block_correct_mean,'m','LineWidth',2);
    %title(sprintf('Average Distances all correct trials across BCI vs. Block Sham: %s (session = %d)',monkey, s));
    xlabel('Distance from Target','FontSize',fontsize)
    ylabel('# of Trials','FontSize',fontsize)
    xlim([0.8 1.93])
    ylim([0 210])
    legend('BCI','Block Sham',sprintf('Mean = %.3f',bci_correct_mean), sprintf('Mean = %.3f',block_correct_mean));
    saveas(e,sprintf('hist_figs/%s/bci_block-correct/hist-correct-mean-s%d.jpg',monkey,s));
    
    %% iterate and concatenate all the vectors for correct vs error into a
    %single histogram and plot the others
    
    fontsize = 16;
    
    if not(isfolder(sprintf('hist_figs/%s/all-correct_error',monkey)))
        mkdir(sprintf('hist_figs/%s/all-correct_error',monkey))
    end
    
    h=figure('Visible','off');
    hold on
    hh1 = histogram(bci_correct_all_mean,'BinWidth',binwidth,'FaceColor',blue);
    hh2 = histogram(bci_error_all_mean,'BinWidth',binwidth,'FaceColor',orange);
    xline(bci_correct_mean,'Color',blue,'LineWidth',2);
    xline(bci_error_mean,'Color',orange,'LineWidth',2);
    %title(sprintf('BCI: Correct vs. Missed Distances Averaged For Each Trial: %s (session = %d)', monkey,s));
    xlabel('Distance from Target','FontSize',fontsize)
    ylabel('# of Trials','FontSize',fontsize)
    xlim([0.8 2])
    ylim([0 210])
    legend('Correct','Missed',sprintf('Mean = %.3f',bci_correct_mean), sprintf('Mean = %.3f',bci_error_mean));
    saveas(h,sprintf('hist_figs/%s/all-correct_error/hist-bci-mean-s%d.jpg',monkey,s));

    p=figure('Visible','off');
    hold on
    histogram(sham_correct_all_mean,'BinWidth',binwidth,'FaceColor',blue);
    histogram(sham_error_all_mean,'BinWidth',binwidth,'FaceColor',orange);
    xline(sham_correct_mean,'Color',blue,'LineWidth',2);
    xline(sham_error_mean,'Color',orange,'LineWidth',2);
    %title(sprintf('BCI Sham: Correct vs. Missed Distances Averaged For Each Trial: %s (session = %d)',monkey,s));
    xlabel('Distance from Target','FontSize',fontsize)
    ylabel('# of Trials','FontSize',fontsize)
    xlim([0.8 2])
    ylim([0 210])
    legend('Correct','Missed',sprintf('Mean = %.3f',sham_correct_mean), sprintf('Mean = %.3f',sham_error_mean));
    saveas(p,sprintf('hist_figs/%s/all-correct_error/hist-sham-mean-s%d.jpg',monkey,s));

    q=figure('Visible','off');
    hold on
    histogram(block_correct_all_mean,'BinWidth',binwidth,'FaceColor',blue);
    histogram(block_error_all_mean,'BinWidth',binwidth,'FaceColor',orange);
    xline(block_correct_mean,'Color',blue,'LineWidth',2);
    xline(block_error_mean,'Color',orange,'LineWidth',2);
    %title(sprintf('Block Sham: Correct vs. Missed Distances Averaged for Each Trial: %s (session = %d)',monkey,s));
    xlabel('Distance from Target','FontSize',fontsize)
    ylabel('# of Trials','FontSize',fontsize)
    xlim([0.8 2])
    ylim([0 210])
    legend('Correct','Missed',sprintf('Mean = %.3f',block_correct_mean), sprintf('Mean = %.3f',block_error_mean));
    saveas(q,sprintf('hist_figs/%s/all-correct_error/hist-block-mean-s%d.jpg',monkey,s));

    % compare distances on all trials across bci vs bci sham vs block sham (bci < bci sham < block sham)
%     bci_all_mean = [bci_correct_all_mean , bci_error_all_mean]; 
%     sham_all_mean = [sham_correct_all_mean, sham_error_all_mean]; 
%     block_all_mean = [block_correct_all_mean, block_error_all_mean]; 
% 
%     d = figure('Visible','off');
%     hold on
%     histogram(bci_all_mean,'BinWidth',binwidth,'FaceColor',blue);
%     histogram(block_all_mean,'BinWidth',binwidth,'FaceColor',orange);
%     xlabel('Distance from Target')
%     ylabel('# of Trials')
%     title(sprintf('Average Distances on all trials across BCI vs. Block Sham (session = %d)',s));
%     legend('BCI','Block Sham')
%     saveas(d,'hist_figs/hist-all-mean-s3.jpg');
end

% Calculating mean across all sessions
bci_correct_all_mean2 = mean(bci_correct_all, 'omitnan');
bci_error_all_mean2 = mean(bci_error_all, 'omitnan');
sham_correct_all_mean2 = mean(sham_correct_all, 'omitnan');
sham_error_all_mean2 = mean(sham_error_all, 'omitnan');
block_correct_all_mean2 = mean(block_correct_all, 'omitnan');
block_error_all_mean2 = mean(block_error_all, 'omitnan');

%% compare distances on error across bci vs block sham

fontsize = 16;

if not(isfolder(sprintf('hist_figs/%s/bci_block-error',monkey)))
    mkdir(sprintf('hist_figs/%s/bci_block-error',monkey))
end

% parameters & colors for plots
binwidth = 0.08;
blue = [0, 0.4470, 0.7410];
orange = [0.8500 0.3250 0.0980];
green = [0.4660, 0.6740, 0.1880];	
gray = [0.17 0.17 0.17];

f = figure('Visible','off');
hold on
h1 = histogram(bci_error_all,'BinWidth',binwidth,'FaceColor','k');
h2 = histogram(block_error_all,'BinWidth',binwidth,'FaceColor','m');
xline(bci_error_all_mean2,'k','LineWidth',2);
xline(block_error_all_mean2,'m','LineWidth',2);
%title(sprintf('Average Distances all missed trials across sessions for BCI vs. Block Sham: %s',monkey));
xlabel('Distance from Target','FontSize',fontsize)
ylabel('# of Trials','FontSize',fontsize)
xlim([0.5 2.5])
ylim([0 2500])
legend('BCI','Block Sham',sprintf('Mean = %.3f',bci_error_all_mean2), sprintf('Mean = %.3f',block_error_all_mean2));
saveas(f,sprintf('hist_figs/%s/bci_block-error/hist-error-mean-all.jpg',monkey));

%% compare distances on correct across bci vs block sham
if not(isfolder(sprintf('hist_figs/%s/bci_block-correct',monkey)))
    mkdir(sprintf('hist_figs/%s/bci_block-correct',monkey))
end

e = figure('Visible','off');
hold on
h3 = histogram(bci_correct_all,'BinWidth',binwidth,'FaceColor','k');
h4 = histogram(block_correct_all,'BinWidth',binwidth,'FaceColor','m');
xline(bci_correct_all_mean2,'k','LineWidth',2);
xline(block_correct_all_mean2,'m','LineWidth',2);
%title(sprintf('Average Distances all correct trials across sessions for BCI vs. Block Sham: %s',monkey));
xlabel('Distance from Target','FontSize',fontsize)
ylabel('# of Trials','FontSize',fontsize)
xlim([0.5 2.5])
ylim([0 2500])
legend('BCI','Block Sham',sprintf('Mean = %.3f',bci_correct_all_mean2), sprintf('Mean = %.3f',block_correct_all_mean2));
saveas(e,sprintf('hist_figs/%s/bci_block-correct/hist-correct-mean-all.jpg',monkey));

%% BCI, BCI Sham, Block Sham all sessions

fontsize = 16;

hall=figure('Visible','off');
hold on
hh1 = histogram(bci_correct_all,'BinWidth',binwidth,'FaceColor',blue);
hh2 = histogram(bci_error_all,'BinWidth',binwidth,'FaceColor',orange);
xline(bci_correct_all_mean2,'Color',blue,'LineWidth',2);
xline(bci_error_all_mean2,'Color',orange,'LineWidth',2);
%title(sprintf('BCI: Correct vs. Missed Distances Averaged For All Sessions: %s', monkey));
xlabel('Distance from Target','FontSize',fontsize)
ylabel('# of Trials','FontSize',fontsize)
xlim([0.5 2.6])
ylim([0 2500])
legend('Correct','Missed',sprintf('Mean = %.3f',bci_correct_all_mean2), sprintf('Mean = %.3f',bci_error_all_mean2));
saveas(hall,sprintf('hist_figs/%s/all-correct_error/hist-bci-mean-all.jpg',monkey));

pall=figure('Visible','off');
hold on
histogram(sham_correct_all,'BinWidth',binwidth,'FaceColor',blue);
histogram(sham_error_all,'BinWidth',binwidth,'FaceColor',orange);
xline(sham_correct_all_mean2,'Color',blue,'LineWidth',2);
xline(sham_error_all_mean2,'Color',orange,'LineWidth',2);
%title(sprintf('BCI Sham: Correct vs. Missed Distances Averaged For All Sessions: %s',monkey));
xlabel('Distance from Target','FontSize',fontsize)
ylabel('# of Trials','FontSize',fontsize)
xlim([0.6 2.4])
ylim([0 300])
legend('Correct','Missed',sprintf('Mean = %.3f',sham_correct_all_mean2), sprintf('Mean = %.3f',sham_error_all_mean2));
saveas(pall,sprintf('hist_figs/%s/all-correct_error/hist-sham-mean-all.jpg',monkey));

qall=figure('Visible','off');
hold on
histogram(block_correct_all,'BinWidth',binwidth,'FaceColor',blue);
histogram(block_error_all,'BinWidth',binwidth,'FaceColor',orange);
xline(block_correct_all_mean2,'Color',blue,'LineWidth',2);
xline(block_error_all_mean2,'Color',orange,'LineWidth',2);
%title(sprintf('Block Sham: Correct vs. Missed Distances Averaged for All Sessions: %s',monkey));
xlabel('Distance from Target','FontSize',fontsize)
ylabel('# of Trials','FontSize',fontsize)
xlim([0.5 2.4])
ylim([0 600])
legend('Correct','Missed',sprintf('Mean = %.3f',block_correct_all_mean2), sprintf('Mean = %.3f',block_error_all_mean2));
saveas(qall,sprintf('hist_figs/%s/all-correct_error/hist-block-mean-all.jpg',monkey));