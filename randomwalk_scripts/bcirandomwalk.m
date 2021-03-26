clear all;

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
    
end

%% Estimating deltaX = eps_t from data to get distribution
bci_c_mean = mean(bci_correct_mean_tot,'omitnan');
bci_e_mean = mean(bci_error_mean_tot,'omitnan');
block_c_mean = mean(block_correct_mean_tot,'omitnan');
block_e_mean = mean(block_error_mean_tot,'omitnan');

bci_c_eps = diff(bci_correct_mean_tot,1,2);
bci_c_eps_mu = mean(mean(bci_c_eps,'omitnan'));
bci_c_eps_std = mean(std(bci_c_eps,'omitnan'));
bci_e_eps = diff(bci_error_mean_tot,1,2);
bci_e_eps_mu = mean(mean(bci_e_eps, 'omitnan'));
bci_e_eps_std = mean(std(bci_e_eps));

block_c_eps = diff(block_correct_mean_tot,1,2);
block_c_eps_mu = mean(mean(block_c_eps,'omitnan'));
block_c_eps_std = mean(std(block_c_eps,'omitnan'));
block_e_eps = diff(block_error_mean_tot,1,2);
block_e_eps_mu = mean(mean(block_e_eps,'omitnan'));
block_e_eps_std = mean(std(block_e_eps));

numberOfSteps = 75; %timesteps

%% BCI Correct
mu_bci_c = bci_c_eps_mu-0.0002; %bias
sigma_bci_c = bci_c_eps_std-0.055; %noisiness
y0_bci_c = bci_c_mean(1); %starting point

% Step sizes
deltax = ones(1, numberOfSteps);
deltay = normrnd(mu_bci_c, sigma_bci_c,[1, numberOfSteps]);
% Starting location
xx = zeros(1, numberOfSteps-1);
x = [1 xx];
yy = zeros(1, numberOfSteps-1);
y = [y0_bci_c yy];

a=figure;
hold on;
% Stepping/walking in random directions
for step = 2 : numberOfSteps
	x(step) = x(step-1) + deltax(step);
	y(step) = y(step-1) + deltay(step);
	plot(x(1:step), y(1:step), 'k-');
end

title(sprintf('Average Distance Over Correct Trials Across Sessions: %s',monkey));
xlim([0 76])
xlabel('Time from BCI start (1 bin = 50 ms)');
ylabel('Avg. Distance');

% Block Correct
mu_block_c = block_c_eps_mu+0.0015; %bias
sigma_block_c = block_c_eps_std-0.124; %noisiness
y0_block_c = block_c_mean(1); %starting point

% Step sizes
deltax = ones(1, numberOfSteps);
deltay = normrnd(mu_block_c, sigma_block_c,[1, numberOfSteps]);
% Starting location
xx = zeros(1, numberOfSteps-1);
x = [1 xx];
yy = zeros(1, numberOfSteps-1);
y = [y0_block_c yy];

hold on;
% Stepping/walking in random directions
for step = 2 : numberOfSteps
	x(step) = x(step-1) + deltax(step);
	y(step) = y(step-1) + deltay(step);
	plot(x(1:step), y(1:step), 'm-');
end

h(1) = plot(1, y0_bci_c,'k-'); hold on;
h(2) = plot(1, y0_block_c,'m-'); 

legend(h, 'BCI','Block Sham','Location','Best');
xlim([0 76])
ylim([0.8 1.4])
saveas(a,sprintf('randomwalk/%s/fake/correct-all.png',monkey)); 
hold off;

%% BCI Missed
mu_bci_e = bci_e_eps_mu+0.00047; %bias
sigma_bci_e = bci_e_eps_std-0.002; %noisiness
y0_bci_e = bci_e_mean(1); %starting point

% Step sizes
deltax = ones(1, numberOfSteps);
deltay_all = [];
for i = 1:21
    y_vals = normrnd(mu_bci_e, sigma_bci_e,[1, numberOfSteps]);
    deltay_all = [deltay_all; y_vals];
end
deltay = mean(deltay_all);
% Starting location
xx = zeros(1, numberOfSteps-1);
x = [1 xx];
yy = zeros(1, numberOfSteps-1);
y = [y0_bci_e yy];

b=figure;
hold on;
% Stepping/walking in random directions
for step = 2 : numberOfSteps
	x(step) = x(step-1) + deltax(step);
	y(step) = y(step-1) + deltay(step);
	plot(x(1:step), y(1:step), 'k-');
end

% Block Missed
mu_block_e = block_e_eps_mu+0.00214; %bias
sigma_block_e = block_e_eps_std-0.004; %noisiness
y0_block_e = block_e_mean(1); %starting point

% Step sizes
deltax = ones(1, numberOfSteps);
deltay_all = [];
for i = 1:21
    y_vals = normrnd(mu_block_e, sigma_block_e,[1, numberOfSteps]);
    deltay_all = [deltay_all; y_vals];
end
deltay = mean(deltay_all);
% Starting location
xx = zeros(1, numberOfSteps-1);
x = [1 xx];
yy = zeros(1, numberOfSteps-1);
y = [y0_block_e yy];

hold on;
% Stepping/walking in random directions
for step = 2 : numberOfSteps
	x(step) = x(step-1) + deltax(step);
	y(step) = y(step-1) + deltay(step);
	plot(x(1:step), y(1:step), 'm-');
end

h(1) = plot(1, y0_bci_e,'k-'); hold on;
h(2) = plot(1, y0_block_e,'m-'); 

title(sprintf('Average Distance Over Missed Trials Across Sessions: %s',monkey));
xlabel('Time from BCI start (1 bin = 50 ms)');
ylabel('Avg. Distance');
legend(h, 'BCI','Block Sham','Location','Best');
xlim([0 76])
ylim([0.8 1.4])
saveas(b,sprintf('randomwalk/%s/fake/error-all.png',monkey)); 