clear;
clear figs;

monkey = 'Pepe';
%monkey = 'Satchel';

load(sprintf('%s/bcidistance.mat',monkey));
load(sprintf('%s/bciacqtime.mat',monkey)); 
load(sprintf('%s/blockshamacqtime.mat',monkey));
load(sprintf('%s/blockshamdistance.mat',monkey));
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
    block_sem_e_tot = [block_sem_e_tot; block_e_sem];
    
end

bci_c_mean = mean(bci_correct_mean_tot,'omitnan');
bci_e_mean = mean(bci_error_mean_tot,'omitnan');
block_c_mean = mean(block_correct_mean_tot,'omitnan');
block_e_mean = mean(block_error_mean_tot,'omitnan');

bci_c_eps = diff(bci_c_mean);
bci_c_eps_mu = mean(bci_c_eps);
bci_c_eps_std = std(bci_c_eps);
bci_e_eps = diff(bci_e_mean);
bci_e_eps_mu = mean(bci_e_eps);
bci_e_eps_std = std(bci_e_eps);

block_c_eps = diff(block_c_mean);
block_c_eps_mu = mean(block_c_eps);
block_c_eps_std = std(block_c_eps);
block_e_eps = diff(block_e_mean);
block_e_eps_mu = mean(block_e_eps);
block_e_eps_std = std(block_e_eps);
