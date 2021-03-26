clear;
clear figs;

%What I think I have made:
%Empirical CDF of acquisition times for correct trials

%Loading in acquisition times
bci = load('bciacqtime.mat'); %148x90
bcisham = load('bcishamacqtime'); %148x10
blocksham = load('blockshamacqtime'); %148x20

%Calibration data? 

%Grabbing the first column (session #1)
bci_s1 = bci.bciacqtime(:,1);
bcisham_s1 = bcisham.bcishamacqtime(:,1);
blocksham_s1 = blocksham.blockshamacqtime(:,1);

%Tossing out the NaN (mistrial) values
bci_s1_c = (bci_s1(~isnan(bci_s1)));
bcisham_s1_c = (bcisham_s1(~isnan(bcisham_s1)));
blocksham_s1_c = (blocksham_s1(~isnan(blocksham_s1)));

%CDF plots
figure
hold on
[f1_bci,x1] = ecdf(bci_s1_c);
[f1_bcisham,x2] = ecdf(bcisham_s1_c);
[f1_blocksham,x3] = ecdf(blocksham_s1_c);
plot(x1,f1_bci,'r')
plot(x2,f1_bcisham, 'b')
plot(x3,f1_blocksham, 'm')
title('Empirical CDF of Acquisition Times')
xlabel('X')
ylabel('F(X)')
legend('BCI','BCI Sham','Block Sham','Location','Best')
