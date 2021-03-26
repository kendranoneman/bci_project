clear all;

monkey = 'Pepe';
%monkey = 'Satchel';

load(sprintf('%s/bcidistance.mat',monkey));
load(sprintf('%s/bciacqtime.mat',monkey)); 
load(sprintf('%s/bcishamdistance.mat',monkey));
load(sprintf('%s/bcishamacqtime.mat',monkey));
load(sprintf('%s/blockshamacqtime.mat',monkey));
load(sprintf('%s/blockshamdistance.mat',monkey));
load(sprintf('%s/sessindindex.mat',monkey));

unique_sessidx = unique(sessindindex);

bci_c_count = []; bci_e_count = [];
sham_c_count = []; sham_e_count = [];
block_c_count = []; block_e_count = [];

for i = 1:length(unique_sessidx)
    s = unique_sessidx(i);
    
    bcidists = bcidistance{i}(:,1:80);
    shamdists = bcishamdistance{i};
    blockdists = blockshamdistance{i};
    
    bci_acqtimes = bciacqtime(sessindindex==s,:); 
    bci_acqtimes = bci_acqtimes(:,1:80); 
    sham_acqtimes = bcishamacqtime(sessindindex==s,:);
    block_acqtimes = blockshamacqtime(sessindindex==s,:); 

    bci_correct_idx = ~isnan(bci_acqtimes);
    sham_correct_idx = ~isnan(sham_acqtimes);
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
    
    sham_trim = [];
    for q3 = 1:size(shamdists,1)
        for m3 = 1:size(shamdists,2)
            aa3 = sham_acqtimes(q3,m3);
            bb3 = shamdists(q3,m3);
            if isnan(aa3)
                sham_trim_dist2=nan;
            else
                mat3 = bb3{1}([1:aa3]);
                sham_trim_dist2 = mat2cell(mat3,1,aa3);
                sham_trim = [sham_trim; sham_trim_dist2];
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
    bci_c_count = [bci_c_count, size(bci_correct_dists,1)];
    bci_error_dists = bcidists(~bci_correct_idx); 
    bci_e_count = [bci_e_count, size(bci_error_dists,1)];
    
    sham_correct_dists = sham_trim;
    sham_c_count = [sham_c_count, size(sham_correct_dists,1)];
    sham_error_dists = shamdists(~sham_correct_idx);
    sham_e_count = [sham_e_count, size(sham_error_dists,1)];
    
    block_correct_dists = block_trim;
    block_c_count = [block_c_count, size(block_correct_dists,1)];
    block_error_dists = blockdists(~block_correct_idx);
    block_e_count = [block_e_count, size(block_error_dists,1)];
end

bci_per_correct = (bci_c_count)./(bci_c_count + bci_e_count) * 100;
bci_per_avg = mean(bci_per_correct);
sham_per_correct = (sham_c_count)./(sham_c_count + sham_e_count) * 100;
sham_per_avg = mean(sham_per_correct);
block_per_correct = (block_c_count)./(block_c_count + block_e_count) * 100;
block_per_avg = mean(block_per_correct);

%% Number of correct/missed in each session
a=figure;
hold on;
x = 1:21;
scatter(x,bci_c_count,'k','filled');
scatter(x,bci_e_count,'k');
scatter(x,sham_c_count, 'b', 'filled');
scatter(x,sham_e_count, 'b');
scatter(x,block_c_count, 'm', 'filled');
scatter(x,block_e_count, 'm');
xlim([1 21]);
xticks(1:21);
xlabel('Session #')
ylabel('Number of Sessions')
legend('BCI Correct','BCI Missed','Sham Correct','Sham Missed', 'Block Correct','Block Missed');
title('Number of Correct & Missed Trials per Session');
saveas(a,sprintf('per_correct/%s-cm-trials.png',monkey));

%% Percent Correct in each session
b=figure;
hold on;
x = 1:21;
scatter(x,bci_per_correct,'k','filled');
scatter(x,sham_per_correct, 'b', 'filled');
scatter(x,block_per_correct, 'm', 'filled');
yline(bci_per_avg, 'k');
yline(sham_per_avg, 'b');
yline(block_per_avg, 'm');
xlim([1 21]);
xticks(1:21);
xlabel('Session #')
ylabel('% Correct Trials')
legend(sprintf('BCI = %.3f',bci_per_avg),sprintf('Sham = %.3f', sham_per_avg),sprintf('Block = %.3f', block_per_avg));
title('Percent Correct per Session');
saveas(b,sprintf('per_correct/%s-per-trials.png',monkey));