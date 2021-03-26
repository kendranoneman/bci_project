clear;
clear figs;

%What I think I have made:
%Empirical CDF of acquisition times for correct trials
%ECDF is defined as the proportion of X values less than or equal to x

%Loading in acquisition times
bci = load('bciacqtime.mat'); %148x90
bcisham = load('bcishamacqtime'); %148x10
blocksham = load('blockshamacqtime'); %148x20
%Calibration data? 

%Grabbing data from a session
s = 10;
bci_s = bci.bciacqtime(:,s);
bcisham_s = bcisham.bcishamacqtime(:,s);
blocksham_s = blocksham.blockshamacqtime(:,s);

%Tossing out the NaN (mistrial) values
bci_s_c = (bci_s(~isnan(bci_s)));
bcisham_s_c = (bcisham_s(~isnan(bcisham_s)));
blocksham_s_c = (blocksham_s(~isnan(blocksham_s)));

%CDF plots (individual trials)
figure
hold on 
[f1_bci,x1] = ecdf(bci_s_c);
[f1_bcisham,x2] = ecdf(bcisham_s_c);
[f1_blocksham,x3] = ecdf(blocksham_s_c);
plot(x1,f1_bci,'r')
plot(x2,f1_bcisham, 'b')
plot(x3,f1_blocksham, 'm')
title('Empirical CDF of Acquisition Times')
xlabel('X')
ylabel('F(X)')
legend('BCI','BCI Sham','Block Sham','Location','Best')