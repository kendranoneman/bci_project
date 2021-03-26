clear;
clear figs;

monkey = 'Pepe';
%monkey = 'Satchel';

load(sprintf('%s/bcidistance.mat',monkey));
load(sprintf('%s/bciacqtime.mat',monkey)); 
load(sprintf('%s/blockshamacqtime.mat',monkey));
load(sprintf('%s/blockshamdistance.mat',monkey));
load(sprintf('%s/sessindindex.mat',monkey));
load(sprintf('%s/thresholds.txt',monkey));

unique_sessidx = unique(sessindindex);
if not(isfolder('psth'))
    mkdir('psth')
end

unique_sessidx = unique(sessindindex);

bci_block_c_mean_diff_tot = [];
bci_block_e_mean_diff_tot = [];
bci_block_c_sem_diff_tot = [];
bci_block_e_sem_diff_tot = [];

bci_sem_c_tot = [];
bci_sem_e_tot = [];
block_sem_c_tot = [];
block_sem_e_tot = [];

bci_correct_mean_tot = [];
bci_error_mean_tot = [];
block_correct_mean_tot = [];
block_error_mean_tot = [];

bci_count_tot = [];
bci_count2_tot = [];
block_count_tot = [];
block_count2_tot = [];

for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);
    t = thresholds(i);
    
    bcidists = bcidistance{i};
    bcidists = bcidists(:,1:80);
    blockdists = blockshamdistance{i}; 
    
    bci_acqtimes = bciacqtime(sessindindex==s,:); 
    bci_acqtimes = bci_acqtimes(:,1:80); 
    block_acqtimes = blockshamacqtime(sessindindex==s,:); 

    bci_correct_idx = ~isnan(bci_acqtimes);
    block_correct_idx = ~isnan(block_acqtimes);
   
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
    block_correct_dists = block_trim;
    block_error_dists = blockdists(~block_correct_idx);
    block_error_dists = block_error_dists(~cellfun('isempty',block_error_dists));
     
    bci_correct_all = [];
    for j = 1:length(bci_correct_dists)
       bci_size = size(bci_correct_dists{j},2);
       num_nan1 = 75-bci_size;
       nan_vect1 = NaN(1,num_nan1);
       bci_nan = [bci_correct_dists{j}, nan_vect1];
       bci_correct_all = [bci_correct_all; bci_nan];
    end
    bci_correct_mean = mean(bci_correct_all,'omitnan');
    bci_correct_mean_tot = [bci_correct_mean_tot; bci_correct_mean];
    nan1 = ~isnan(bci_correct_all);
    bci_count = sum(nan1);
    bci_count_tot = [bci_count_tot; bci_count];
    bci_c_sem = std(bci_correct_all,'omitnan')./sqrt(bci_count);
    bci_sem_c_tot = [bci_sem_c_tot; bci_c_sem];
    
    block_correct_all = [];
    for j = 1:length(block_correct_dists)
       block_size = size(block_correct_dists{j},2);
       num_nan2 = 75-block_size;
       nan_vect2 = NaN(1,num_nan2);
       block_nan = [block_correct_dists{j}, nan_vect2];
       block_correct_all = [block_correct_all; block_nan];
    end
    block_correct_mean = mean(block_correct_all,'omitnan');
    block_correct_mean_tot = [block_correct_mean_tot; block_correct_mean];
    nan2 = ~isnan(block_correct_all);
    block_count = sum(nan2);
    block_count_tot = [block_count_tot; block_count];
    block_c_sem = std(block_correct_all,'omitnan')./sqrt(block_count);
    block_sem_c_tot = [block_sem_c_tot; block_c_sem];
    
    bci_block_c_mean_diff_tot = [bci_block_c_mean_diff_tot ; bci_correct_mean - block_correct_mean];
    bci_block_c_sem_diff_tot = [bci_block_c_sem_diff_tot ; bci_c_sem - block_c_sem];
    
    bci_error_all = [];
    for j = 1:length(bci_error_dists)
       bci_size = size(bci_error_dists{j},2);
       num_nan1 = 75-bci_size;
       nan_vect1 = NaN(1,num_nan1);
       bci_nan = [bci_error_dists{j}, nan_vect1];
       bci_error_all = [bci_error_all; bci_nan];
    end
    bci_error_mean = mean(bci_error_all,'omitnan');
    bci_error_mean_tot = [bci_error_mean_tot; bci_error_mean];
    nan3 = ~isnan(bci_error_all);
    bci_count2 = sum(nan3);
    bci_count2_tot = [bci_count2_tot; bci_count2];
    bci_e_sem = std(bci_error_all,'omitnan')./sqrt(bci_count2);
    bci_sem_e_tot = [bci_sem_e_tot; bci_e_sem];
    
    block_error_all = [];
    trial_avg = [];
    for j = 1:length(block_error_dists)
       block_size = size(block_error_dists{j},2);
       num_nan2 = 75-block_size;
       nan_vect2 = NaN(1,num_nan2);
       block_nan = [block_error_dists{j}, nan_vect2];
       block_error_all = [block_error_all; block_nan];
       trial_avg = [trial_avg; mean(block_error_dists{j})];
    end
    block_error_mean = mean(block_error_all,'omitnan');
    block_error_mean_tot = [block_error_mean_tot; block_error_mean];
    nan4 = ~isnan(block_error_all);
    block_count2 = sum(nan4);
    block_count2_tot = [block_count2_tot; block_count2];
    block_e_sem = std(block_error_all,'omitnan')./sqrt(block_count2);
    block_sem_e_tot = [block_sem_e_tot; block_e_sem];
    
    bci_block_e_mean_diff_tot = [bci_block_e_mean_diff_tot ; bci_error_mean - block_error_mean];
    bci_block_e_sem_diff_tot = [bci_block_e_sem_diff_tot ; bci_e_sem - block_e_sem];
    
%% Histogram of raw data for each trial
    
    binwidth = 0.02;
    blue = [0, 0.4470, 0.7410];
    orange = [0.8500 0.3250 0.0980];
    green = [0.4660, 0.6740, 0.1880];	
    gray = [0.17 0.17 0.17];

%     max_1 = max(bci_correct_mean); min_1 = min(bci_correct_mean);
%     max_2 = max(bci_error_mean); min_2 = min(bci_error_mean);
%     max_3 = max(block_correct_mean); min_3 = min(block_correct_mean);
%     max_4 = max(block_error_mean); min_4 = min(block_error_mean);
%     
%     maxes = [max_1, max_2, max_3, max_4];
%     mins = [min_1, min_2, min_3, min_4];
%     maxmax = max(maxes);
%     minmin = min(mins);
%     
%     a = figure('Visible','Off'); 
%     subplot(2,2,1)
%     histogram(bci_correct_mean,'BinWidth',binwidth,'FaceColor',blue);
%     title(sprintf('BCI Correct (s = %d): %s',s,monkey));
%     xlabel('Distances')
%     ylabel('Number of Time Bins')
%     xlim([minmin-0.1 maxmax+0.1])
%     ylim([0 30])
%     
%     subplot(2,2,2)
%     histogram(bci_error_mean,'BinWidth',binwidth,'FaceColor',orange);
%     title(sprintf('BCI Missed (s = %d): %s',s,monkey));
%     xlabel('Distances')
%     ylabel('Number of Time Bins')
%     xlim([minmin-0.1 maxmax+0.1])
%     ylim([0 30])
%     
%     subplot(2,2,3)
%     histogram(block_correct_mean,'BinWidth',binwidth,'FaceColor',blue);
%     title(sprintf('Block Correct (s = %d): %s',s,monkey));
%     xlabel('Distances')
%     ylabel('Number of Time Bins')
%     xlim([minmin-0.1 maxmax+0.1])
%     ylim([0 30])
%     
%     subplot(2,2,4)
%     histogram(block_error_mean,'BinWidth',binwidth,'FaceColor',orange);
%     title(sprintf('Block Missed (s = %d): %s',s,monkey));
%     xlabel('Distances')
%     ylabel('Number of Time Bins')
%     xlim([minmin-0.1 maxmax+0.1])
%     ylim([0 30])
%     
%     saveas(a,sprintf('hists/%s/hist-s%d.png',monkey,s));
    
    %% Plot Block Sham Missed Trials PSTH
%     x = (1:75)*50;
%     for j = 1:length(block_error_dists)
%        b=figure('Visible','Off');
%        plot(x,block_error_dists{j})
%        title(sprintf('Distance for Missed Trials vs. Time: %s (s=%d, t=%d)',monkey, s,j));
%        xlabel('Time from BCI start (ms)');
%        ylabel('Distance');
%        saveas(b,sprintf('hists/%s/block/trials/solo/psth-s%d-t%d.png',monkey,s,j));
%     end
%     
%     for jj = 1:length(block_error_dists)
%         c = figure('Visible','Off');
%         histogram(block_error_dists{jj})
%         title(sprintf('Distance for Missed Trials vs. Time: %s (s=%d, t=%d)',monkey,s,jj));
%         xlabel('Distance')
%         saveas(c,sprintf('hists/%s/block/trials/solo/hist-s%d-t%d.png',monkey,s,jj));
%     end

    %% Plot avg. distance for each trial in Block Sham Missed
       xx = 1:size(block_error_dists,1);
       dd = figure('Visible','Off');
       scatter(xx, trial_avg,'filled');
       [M3, I3] = maxk(trial_avg, 3)
       dx = 1.5; dy = 0.1;
       text(I3(1)+dx,M3(1),sprintf('t = %d',I3(1))); 
       text(I3(2)+dx,M3(2),sprintf('t = %d',I3(2)));
       text(I3(3)+dx,M3(3),sprintf('t = %d',I3(3))); 
       title(sprintf('Distance for Missed Block Sham Trials vs. Time: %s (s=%d)',monkey,s));
       xlabel('Trial Number');
       xticks(1:10:size(block_error_dists,1))
       ylim([0.6 2.5])
       ylabel('Average Distance in Trial');
       saveas(dd,sprintf('hists/%s/block/trials/solo/avg-s%d.png',monkey,s));
       
       binwidth = 0.02
       for ii = 1:3
           I = I3(ii)
           dd2 = figure('Visible','Off');
           x = (1:75)*50;
           plot(x,block_error_dists{I})
           title(sprintf('Distance for Missed Block Sham Trials vs. Time: %s (s=%d, t=%d)',monkey, s,I));
           xlabel('Time from BCI start (ms)');
           ylabel('Distance');
           ylim([0 12])
           saveas(dd2,sprintf('hists/%s/block/trials/solo/giveup/psth-s%d-t%d.png',monkey,s,I)); 
           
           dd3 = figure('Visible','Off');
           histogram(block_error_dists{I},'BinWidth',binwidth)
           title(sprintf('Distance for Missed Block Sham Trials: %s (s=%d, t=%d)',monkey, s,I));
           xlabel('Distances');
           ylabel('Number of Time Bins');
           %ylim([0 12])
           saveas(dd2,sprintf('hists/%s/block/trials/solo/giveup/hist-s%d-t%d.png',monkey,s,I));
       end
end

%% All sessions
bci_c_mean = mean(bci_correct_mean_tot,'omitnan');
bci_e_mean = mean(bci_error_mean_tot,'omitnan');
block_c_mean = mean(block_correct_mean_tot,'omitnan');
block_e_mean = mean(block_error_mean_tot,'omitnan');

bci_count_mean = mean(bci_count_tot);
bci_count2_mean = mean(bci_count2_tot);
block_count_mean = mean(block_count_tot);
block_count2_mean = mean(block_count2_tot);

bci_c_sem_mean = std(bci_correct_mean_tot,'omitnan')./sqrt(bci_count_mean);
bci_e_sem_mean = std(bci_error_mean_tot,'omitnan')./sqrt(bci_count2_mean);
block_c_sem_mean = std(bci_error_mean_tot,'omitnan')./sqrt(block_count_mean);
block_e_sem_mean = std(block_error_mean_tot,'omitnan')./sqrt(block_count2_mean);

bci_block_c_mean_diff = mean(bci_block_c_mean_diff_tot);
bci_block_c_sem_diff = mean(bci_block_c_sem_diff_tot);
bci_block_e_mean_diff = mean(bci_block_e_mean_diff_tot);
bci_block_e_sem_diff = mean(bci_block_e_sem_diff_tot);


%% Histogram of raw data for all sessions
% h1 = histfit(bci_c_mean',10);
% h1(1).FaceColor = [0.6 0.6 0.6];
% h1(2).Color = [0.2 0.2 0.2];
% pd = fitdist(bci_c_mean','Normal')
% xline(pd.mu,'k--', 'LineWidth', 2);
% hold on;
% 
% h2 = histfit(block_c_mean',10);
% h2(1).FaceColor = [251,227,251]./255;
% h2(2).Color = [1, 0, 1];
% pd2 = fitdist(block_c_mean', 'Normal')
% xline(pd2.mu,'m--', 'LineWidth', 2);