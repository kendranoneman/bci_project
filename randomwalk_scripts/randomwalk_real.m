clear; clear figs;

monkey = 'Pepe'; 
%monkey = 'Satchel';

load(sprintf('%s/thresholds.txt',monkey));
load(sprintf('%s/sessindindex.mat',monkey));

load(sprintf('%s/bcidistance.mat',monkey));
load(sprintf('%s/bciacqtime.mat',monkey)); 
load(sprintf('%s/blockshamacqtime.mat',monkey));
load(sprintf('%s/blockshamdistance.mat',monkey));
load(sprintf('%s/bcishamacqtime.mat',monkey));
load(sprintf('%s/bcishamdistance.mat',monkey));

unique_sessidx = unique(sessindindex);

tot_avg_dist = [];

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
    %t = thresholds(i);
    
    bcidists = bcidistance{i}(:,1:80);
    shamdists = bcishamdistance{i}; 
    shamall = shamdists(:); shamall2 = [];
    for j1 = 1:length(shamall)
        d1 = shamall{j1};
        shamall2 = [shamall2, d1];
    end
    blockdists = blockshamdistance{i}; 
    blockall = blockdists(:); blockall2 = [];
    for k1 = 1:length(blockall)
        f1 = blockall{k1};
        blockall2 = [blockall2, f1];
    end
    dists = [shamall2, blockall2];
    avg_dist = mean(dists);
    tot_avg_dist = [tot_avg_dist, avg_dist];
    
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
       bci_nan = [bci_correct_dists{j}, nan_vect1];
       bci_correct_all = [bci_correct_all; bci_nan];
    end
    bci_correct_mean = mean(bci_correct_all,'omitnan');
    bci_correct_mean_tot = [bci_correct_mean_tot; bci_correct_mean];
    nan1 = ~isnan(bci_correct_all);
    bci_count = sum(nan1);
    bci_count_tot = [bci_count_tot; bci_count];
    bci_c_sem = std(bci_correct_all,'omitnan')./sqrt(bci_count);
    
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
    
    %bci_block_c_mean_diff_tot = [bci_block_c_mean_diff_tot ; bci_correct_mean - block_correct_mean];
    %bci_block_c_sem_diff_tot = [bci_block_c_sem_diff_tot ; bci_c_sem - block_c_sem];
    
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
    
    block_error_all = [];
    for j = 1:length(block_error_dists)
       block_size = size(block_error_dists{j},2);
       num_nan2 = 75-block_size;
       nan_vect2 = NaN(1,num_nan2);
       block_nan = [block_error_dists{j}, nan_vect2];
       block_error_all = [block_error_all; block_nan];
    end
    block_error_mean = mean(block_error_all,'omitnan');
    block_error_mean_tot = [block_error_mean_tot; block_error_mean];
    nan4 = ~isnan(block_error_all);
    block_count2 = sum(nan4);
    block_count2_tot = [block_count2_tot; block_count2];
    block_e_sem = std(block_error_all,'omitnan')./sqrt(block_count2);

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

%% Plotting PSTH for Correct Across All Sessions
t=figure('Visible','Off');
blue = [0, 0.4470, 0.7410];
orange = [0.8500 0.3250 0.0980];
hold on 
x = (1:75);
hold on
%diff_yu1 = bci_c_mean+bci_c_sem_mean; diff_yu1(isnan(diff_yu1))=0;
%diff_yl1 = bci_c_mean-bci_c_sem_mean; diff_yl1(isnan(diff_yl1))=0;
%diff_yu2 = block_c_mean+block_c_sem_mean; diff_yu2(isnan(diff_yu2))=0;
%diff_yl2 = block_c_mean-block_c_sem_mean; diff_yl2(isnan(diff_yl2))=0;
%fill([x fliplr(x)], [diff_yu1 fliplr(diff_yl1)], [.9 .9 .9], 'linestyle', 'none')
%fill([x fliplr(x)], [diff_yu2 fliplr(diff_yl2)], [251,227,251]./255, 'linestyle','none')
plot(x,bci_c_mean,'k');
plot(x,block_c_mean,'m');
yline(0,'--k')
ylim([0.8 1.4])
xlim([0 76])
title(sprintf('Average Distance Over Correct Trials Across Sessions: %s',monkey));
xlabel('Time from BCI start (1 bin = 50 ms)');
ylabel('Avg. Distance');
legend('BCI','Block Sham','Location','Best');
saveas(t,sprintf('randomwalk/%s/real/correct-all.png',monkey));

%% Plotting PSTH for Error Across All Sessions
t=figure('Visible','Off');
blue = [0, 0.4470, 0.7410];
orange = [0.8500 0.3250 0.0980];
hold on 
x = (1:75);
%diff_yu1 = bci_e_mean+bci_e_sem_mean; diff_yu1(isnan(diff_yu1))=0;
%diff_yl1 = bci_e_mean-bci_e_sem_mean; diff_yl1(isnan(diff_yl1))=0;
%diff_yu2 = block_e_mean+block_e_sem_mean; diff_yu2(isnan(diff_yu2))=0;
%diff_yl2 = block_e_mean-block_e_sem_mean; diff_yl2(isnan(diff_yl2))=0;
%fill([x fliplr(x)], [diff_yu1 fliplr(diff_yl1)], [.9 .9 .9], 'linestyle', 'none')
%fill([x fliplr(x)], [diff_yu2 fliplr(diff_yl2)], [251,227,251]./255, 'linestyle','none')
plot(x,bci_e_mean,'k');
plot(x,block_e_mean,'m');
yline(0,'--k')
ylim([0.8 1.4])
xlim([0 76])
title(sprintf('Average Distance Over Missed Trials Across Sessions: %s',monkey));
xlabel('Time from BCI start (1 bin = 50 ms)');
ylabel('Avg. Distance');
legend('BCI','Block Sham','Location','Best');
saveas(t,sprintf('randomwalk/%s/real/error-all.png',monkey));    
