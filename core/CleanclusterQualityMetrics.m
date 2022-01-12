
clear
cd('F:\Data\sfm\QualityMeasures\0095');                                     %Enter dir to analyze
load('clusterQualityMetrics.mat')

analytic_matrix = [];
analytic_matrix(:, 1) = cids;
analytic_matrix(:, 2) = cgs';
analytic_matrix(:, 3) = cR;
analytic_matrix(:, 4) = isiV;
analytic_matrix(:, 5) = noise_rate;
analytic_matrix(:, 6) = uQ;