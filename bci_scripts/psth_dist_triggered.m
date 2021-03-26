clear;
clear figs;

monkey = 'Pepe';
%monkey = 'Satchel';

load(sprintf('%s/bcidistance.mat',monkey));
load(sprintf('%s/bciacqtime.mat',monkey)); 
load(sprintf('%s/blockshamacqtime.mat',monkey));
load(sprintf('%s/blockshamdistance.mat',monkey));
load(sprintf('%s/bcishamacqtime.mat',monkey));
load(sprintf('%s/bcishamdistance.mat',monkey));
load(sprintf('%s/sessindindex.mat',monkey));

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

bci_count_tot = [];
bci_count2_tot = [];
block_count_tot = [];
block_count2_tot = [];

bci_correct_mean_tot = [];
bci_error_mean_tot = [];
block_correct_mean_tot = [];
block_error_mean_tot = [];

for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);
    
    bcidists = bcidistance{i};
    bcidists = bcidists(:,1:80);
    blockdists = blockshamdistance{i}; 
    
    bci_acqtimes = bciacqtime(sessindindex==s,:); 
    bci_acqtimes = bci_acqtimes(:,1:80); 
    block_acqtimes = blockshamacqtime(sessindindex==s,:); 

    bci_correct_idx = ~isnan(bci_acqtimes);
    block_correct_idx = ~isnan(block_acqtimes);
  
    bci_trim = [];
    for q = 1:size(bcidists,1)
        for m = 1:size(bcidists,2)
            aa2 = bci_acqtimes(q,m);
            bb2 = bcidists(q,m);
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
     
    bci_correct_all = [];
    for j = 1:length(bci_correct_dists)
       bci_size = size(bci_correct_dists{j},2);
       num_nan1 = 75-bci_size;
       nan_vect1 = NaN(1,num_nan1);
       bci_nan = [nan_vect1, bci_correct_dists{j}];
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
       block_nan = [nan_vect2, block_correct_dists{j}];
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
    bci_block_c_sem_diff = std(bci_block_c_mean_diff_tot,'omitnan')./sqrt(bci_count);
    bci_block_c_sem_diff_tot = [bci_block_c_sem_diff_tot ; bci_block_c_sem_diff];
    
    bci_error_all = [];
    for j = 1:length(bci_error_dists)
       bci_size = size(bci_error_dists{j},2);
       num_nan1 = 75-bci_size;
       nan_vect1 = NaN(1,num_nan1);
       bci_nan = [nan_vect1, bci_error_dists{j}];
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
    for j = 1:length(block_error_dists)
       block_size = size(block_error_dists{j},2);
       num_nan2 = 75-block_size;
       nan_vect2 = NaN(1,num_nan2);
       block_nan = [nan_vect2, block_error_dists{j}];
       block_error_all = [block_error_all; block_nan];
    end
    block_error_mean = mean(block_error_all,'omitnan');
    block_error_mean_tot = [block_error_mean_tot; block_error_mean];
    nan4 = ~isnan(block_error_all);
    block_count2 = sum(nan4);
    block_count2_tot = [block_count2_tot; block_count2];
    block_e_sem = std(block_error_all,'omitnan')./sqrt(block_count2);
    block_sem_e_tot = [block_sem_e_tot; block_e_sem];
    
    bci_block_e_mean_diff_tot = [bci_block_e_mean_diff_tot ; bci_error_mean - block_error_mean];
    bci_block_e_sem_diff = std(bci_block_e_mean_diff_tot,'omitnan')./sqrt(bci_count2);
    bci_block_e_sem_diff_tot = [bci_block_e_sem_diff_tot ; bci_block_e_sem_diff];
    
%     %% Print number of trials counted
%     c=figure('Visible','Off');
%     hold on
%     x = (1:75)*50;
%     plot(x,bci_c_sem,'k');
%     plot(x,block_c_sem,'m');
%     title(sprintf('Standard Error vs. Time (session = %d)', s));
%     xlabel('Time from BCI start (ms)');
%     ylabel('Standard Error');
%     legend('BCI','Block Sham','Location','Best');
%     saveas(c,sprintf('psth/psth-c-count-s%d.png',s));    

%% Plotting PSTH for Correct vs. Missed for BCI Trials
%     blue = [0, 0.4470, 0.7410];
%     orange = [0.8500 0.3250 0.0980];
%     a=figure('Visible','Off');
%     hold on 
%     x = (-75:-1)*50;
%     bci_yu = bci_correct_mean+bci_c_sem; bci_yu(isnan(bci_yu))=0;
%     bci_yl = bci_correct_mean-bci_c_sem; bci_yl(isnan(bci_yl))=0;
%     bci2_yu = bci_error_mean+bci_e_sem; bci2_yu(isnan(bci2_yu))=0;
%     bci2_yl = bci_error_mean-bci_e_sem; bci2_yl(isnan(bci2_yl))=0;
%     colororder({'k','k'})
%     yyaxis left
%     fill([x fliplr(x)], [bci_yu fliplr(bci_yl)], [219,233,253]./255, 'linestyle', 'none')
%     fill([x fliplr(x)], [bci2_yu fliplr(bci2_yl)], [253,238,219]./255, 'linestyle','none')
%     plot(x,bci_correct_mean,'-','Color',blue);
%     plot(x,bci_error_mean,'-','Color',orange);
%     %ylim([0.5 1.6])
%     title(sprintf('Average Distance Over BCI Trials vs. Time: %s (session = %d)',monkey, s));
%     xlabel('Time from BCI end (ms)');
%     ylabel('Avg. Distance');
%     yyaxis right
%     plot(x,bci_count,'--','Color',blue);
%     plot(x,bci_count2,'--','Color',orange);
%     %ylim([0 bci_count2(1)+50])
%     ylabel('Number of Trials');
%     legend('Correct','Missed','Location','Best');
%     saveas(a,sprintf('psth/%s/end/psth-trig-cvm-s%d.png',monkey,s));

%% Plotting PSTH for Correct Trials

    fontsize = 16;

    a=figure('Visible','Off');
    hold on 
    x = (-75:-1)*50;
    bci_yu1 = bci_correct_mean+bci_c_sem; bci_yu1(isnan(bci_yu1))=0;
    bci_yl1 = bci_correct_mean-bci_c_sem; bci_yl1(isnan(bci_yl1))=0;
    block_yu1 = block_correct_mean+block_c_sem; block_yu1(isnan(block_yu1))=0;
    block_yl1 = block_correct_mean-block_c_sem; block_yl1(isnan(block_yl1))=0;
    colororder({'k','k'})
    yyaxis left
    fill([x fliplr(x)], [bci_yu1 fliplr(bci_yl1)], [.9 .9 .9], 'linestyle', 'none')
    fill([x fliplr(x)], [block_yu1 fliplr(block_yl1)], [251,227,251]./255, 'linestyle','none')
    plot(x,bci_correct_mean,'k-');
    plot(x,block_correct_mean, 'm-');
    ylim([0.5 3.0])
    %title(sprintf('Average Distance Over Correct Trials vs. Time: %s (session = %d)',monkey, s));
    xlabel('Time from BCI end (ms)','FontSize',fontsize);
    ylabel('Avg. Distance','FontSize',fontsize);
    yyaxis right
    plot(x,bci_count,'k--');
    plot(x,block_count,'--m');
    ylim([0 500])
    ylabel('Number of Trials','FontSize',fontsize);
    legend('BCI','Block Sham','Location','northeast');
    saveas(a,sprintf('psth/%s/end/psth-trig-correct-s%d.png',monkey,s));

     %% Plotting PSTH for Missed Trials
    b=figure('Visible','Off');
    hold on 
    x = (-75:-1)*50;
    bci_yu2 = bci_error_mean+bci_e_sem; bci_yu2(isnan(bci_yu2))=0;
    bci_yl2 = bci_error_mean-bci_e_sem; bci_yl2(isnan(bci_yl2))=0;
    block_yu2 = block_error_mean+block_e_sem; block_yu2(isnan(block_yu2))=0;
    block_yl2 = block_error_mean-block_e_sem; block_yl2(isnan(block_yl2))=0;
    colororder({'k','k'})
    yyaxis left
    fill([x fliplr(x)], [bci_yu2 fliplr(bci_yl2)], [.9 .9 .9], 'linestyle', 'none')
    fill([x fliplr(x)], [block_yu2 fliplr(block_yl2)], [251,227,251]./255, 'linestyle','none')
    plot(x,bci_error_mean,'k-');
    plot(x,block_error_mean, 'm-');
    ylim([0.5 3.0])
    %title(sprintf('Average Distance Over Missed Trials vs. Time: %s (session = %d)',monkey, s));
    xlabel('Time from BCI end (ms)','FontSize',fontsize);
    ylabel('Avg. Distance','FontSize',fontsize);
    yyaxis right
    plot(x,bci_count2,'--k');
    plot(x,block_count2, '--m');
    ylim([0 500])
    ylabel('Number of Trials','FontSize',fontsize);
    legend('BCI','Block Sham','Location','northeast');
    saveas(b,sprintf('psth/%s/end/psth-trig-error-s%d.png',monkey,s));
end
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

bci_count_all = sum(bci_count_tot);
bci_count2_all = sum(bci_count2_tot);
block_count_all = sum(block_count_tot);
block_count2_all = sum(block_count2_tot);

%% Plotting PSTH for Correct vs. Missed Trials
d=figure('Visible','Off');
blue = [0, 0.4470, 0.7410];
orange = [0.8500 0.3250 0.0980];
hold on 
x = (-75:-1)*50;
diff_yu1 = bci_block_c_mean_diff+bci_block_c_sem_diff; diff_yu1(isnan(diff_yu1))=0;
diff_yl1 = bci_block_c_mean_diff-bci_block_c_sem_diff; diff_yl1(isnan(diff_yl1))=0;
diff_yu2 = bci_block_e_mean_diff+bci_block_e_sem_diff; diff_yu2(isnan(diff_yu2))=0;
diff_yl2 = bci_block_e_mean_diff-bci_block_e_sem_diff; diff_yl2(isnan(diff_yl2))=0;
fill([x fliplr(x)], [diff_yu1 fliplr(diff_yl1)], [219,233,253]./255, 'linestyle', 'none')
fill([x fliplr(x)], [diff_yu2 fliplr(diff_yl2)], [253,238,219]./255, 'linestyle','none')
plot(x,bci_block_c_mean_diff,'Color', blue);
plot(x,bci_block_e_mean_diff,'Color',orange);
yline(0,'--k')
ylim([-0.2 0.1])
title(sprintf('Difference in Average Distance Between BCI and Block Sham Over Trials vs. Time: %s',monkey));
xlabel('Time from BCI end (ms)');
ylabel('Difference in Avg. Distance (BCI - Block Sham)');
legend('Correct','Missed','Location','Best');
saveas(d,sprintf('psth/%s/end/psth-trig-diff.png',monkey));

%% Plotting PSTH for Correct Across All Sessions

fontsize = 16;

t=figure('Visible','Off');
blue = [0, 0.4470, 0.7410];
orange = [0.8500 0.3250 0.0980];
hold on 
x = (-75:-1)*50;
diff_yu1a = bci_c_mean+bci_c_sem_mean; diff_yu1a(isnan(diff_yu1a))=0;
diff_yl1a = bci_c_mean-bci_c_sem_mean; diff_yl1a(isnan(diff_yl1a))=0;
diff_yu2a = block_c_mean+block_c_sem_mean; diff_yu2a(isnan(diff_yu2a))=0;
diff_yl2a = block_c_mean-block_c_sem_mean; diff_yl2a(isnan(diff_yl2a))=0;
colororder({'k','k'})
yyaxis left
fill([x fliplr(x)], [diff_yu1a fliplr(diff_yl1a)], [.9 .9 .9], 'linestyle', 'none')
fill([x fliplr(x)], [diff_yu2a fliplr(diff_yl2a)], [251,227,251]./255, 'linestyle','none')
plot(x,bci_c_mean,'k-');
plot(x,block_c_mean,'m-');
if min(block_c_mean) < min(bci_c_mean)
    ymin = min(block_c_mean) - 0.2;
else
    ymin = min(bci_c_mean) - 0.2;
end
if max(block_c_mean) > max(bci_c_mean)
    ymax = max(block_c_mean) + 0.2;
else
    ymax = max(bci_c_mean) + 0.2;
end
ylim([0.3 1.9])
%title(sprintf('Average Distance Over Correct Trials Across Sessions vs. Time: %s',monkey));
xlabel('Time from BCI end (ms)','FontSize',fontsize);
ylabel('Avg. Distance','FontSize',fontsize);
yyaxis right
plot(x,bci_count_all,'k--');
plot(x,block_count_all,'--m');cylim([0 10000])
ylabel('Sum of Trials Across Sessions','FontSize',fontsize);
legend('BCI','Block Sham','Location','Best');
saveas(t,sprintf('psth/%s/end/psth-trig-correct-all.png',monkey));

%% Plotting PSTH for Error Across All Sessions
t=figure('Visible','Off');
blue = [0, 0.4470, 0.7410];
orange = [0.8500 0.3250 0.0980];
hold on 
x = (-75:-1)*50;
diff_yu1b = bci_e_mean+bci_e_sem_mean; diff_yu1b(isnan(diff_yu1b))=0;
diff_yl1b = bci_e_mean-bci_e_sem_mean; diff_yl1b(isnan(diff_yl1b))=0;
diff_yu2b = block_e_mean+block_e_sem_mean; diff_yu2b(isnan(diff_yu2b))=0;
diff_yl2b = block_e_mean-block_e_sem_mean; diff_yl2b(isnan(diff_yl2b))=0;
colororder({'k','k'})
yyaxis left
fill([x fliplr(x)], [diff_yu1b fliplr(diff_yl1b)], [.9 .9 .9], 'linestyle', 'none')
fill([x fliplr(x)], [diff_yu2b fliplr(diff_yl2b)], [251,227,251]./255, 'linestyle','none')
plot(x,bci_e_mean,'k-');
plot(x,block_e_mean,'m-');
ylim([0.3 1.9])
%title(sprintf('Average Distance Over Error Trials Across Sessions vs. Time: %s',monkey));
xlabel('Time from BCI end (ms)','FontSize',fontsize);
ylabel('Avg. Distance','FontSize',fontsize);
yyaxis right
plot(x,bci_count2_all,'k--');
plot(x,block_count2_all,'--m');
ylim([0 10000])
ylabel('Sum of Trials Across Sessions','FontSize',fontsize);
legend('BCI','Block Sham','Location','Best');
saveas(t,sprintf('psth/%s/end/psth-trig-error-all.png',monkey));

%% Plotting PSTH for Correct vs. Missed Across All BCI Sessions
t=figure('Visible','Off');
blue = [0, 0.4470, 0.7410];
orange = [0.8500 0.3250 0.0980];
hold on 
x = (-75:-1)*50;
diff_yu1b = bci_c_mean+bci_c_sem_mean; diff_yu1b(isnan(diff_yu1b))=0;
diff_yl1b = bci_c_mean-bci_c_sem_mean; diff_yl1b(isnan(diff_yl1b))=0;
diff_yu2b = bci_e_mean+bci_e_sem_mean; diff_yu2b(isnan(diff_yu2b))=0;
diff_yl2b = bci_e_mean-bci_e_sem_mean; diff_yl2b(isnan(diff_yl2b))=0;
colororder({'k','k'})
yyaxis left
fill([x fliplr(x)], [diff_yu1b fliplr(diff_yl1b)], [219,233,253]./255, 'linestyle', 'none')
fill([x fliplr(x)], [diff_yu2b fliplr(diff_yl2b)], [253,238,219]./255, 'linestyle','none')
plot(x,bci_c_mean,'-','Color',blue);
plot(x,bci_e_mean,'-','Color',orange);
if min(bci_c_mean) < min(bci_e_mean)
    ymin2 = min(bci_c_mean) - 0.1;
else
    ymin2 = min(bci_e_mean) - 0.1;
end
if max(bci_c_mean) > max(bci_e_mean)
    ymax2 = max(bci_c_mean) + 0.3;
else
    ymax2 = max(bci_e_mean) + 0.3;
end
ylim([0.7 1.6])
title(sprintf('Average Distance Over Correct/Missed Trials Across BCI Sessions vs. Time: %s',monkey));
xlabel('Time from BCI end (ms)');
ylabel('Avg. Distance');
yyaxis right
plot(x,bci_count_all,'--','Color',blue);
plot(x,bci_count2_all,'--','Color',orange);
ylim([0 7000])
ylabel('Sum of Trials Across Sessions');
legend('Correct','Missed','Location','Best');
saveas(t,sprintf('psth/%s/end/psth-trig-cvm-all.png',monkey));

%% Vertical difference b/t BCI & Block
diff_c = block_c_mean - bci_c_mean;
diff_c_avg = mean(diff_c);
diff_e = block_e_mean - bci_e_mean;
diff_e_avg = mean(diff_e);

tt=figure;
hold on;
x = (-75:-1)*50;
plot(x,diff_c,'Color', blue);
plot(x,diff_e, 'Color', orange);
yline(0, '--k');
xlabel('Time from BCI end (ms)');
ylabel('Difference in Avg. Distance (Block - BCI)')
ylim([-0.2 0.3])
legend('Correct', 'Missed');
title(sprintf('Difference in Avg. Distance (Block Sham - BCI) Across Sessions: %s',monkey));
saveas(tt,sprintf('psth/%s/end/psth-diff-lines.png',monkey));

binwidth=0.005;
tt2=figure;
hold on;
histogram(diff_c,'BinWidth',binwidth,'FaceColor',blue);
histogram(diff_e,'BinWidth',binwidth,'FaceColor',orange);
xline(0, '--k','LineWidth',2);
xline(diff_c_avg,'--','Color',blue,'LineWidth',3)
xline(diff_e_avg,'--','Color',orange,'LineWidth',3)
xlim([-0.18 0.3])
ylim([0 16]);
xlabel('Difference in Avg. Distance (Block Sham - BCI)');
ylabel('Number of Time Bins');
legend(sprintf('Correct = %.4f',diff_c_avg), sprintf('Missed = %.4f',diff_e_avg));
title(sprintf('Difference in Avg. Distance (Block Sham - BCI) Across Sessions: %s', monkey));
saveas(tt2,sprintf('psth/%s/end/psth-diff-hist.png',monkey));