function [clusterIDs, unitQuality, contaminationRate, noise_rate] = maskedClusterQualitySparse(clu, fet, fetInds, fetNchans)
% - clu is 1 x nSpikes
% - fet is nSpikes x nPCsPerChan x nInclChans
% - fetInds is nClusters x nInclChans (sorted in descending order of
% relevance for this template)
% - fetN is an integer, the number of features to

if nargin<4
    fetNchans = min(4, size(fetInds,2)); % number of channels to use
end
nFetPerChan = size(fet,2);
fetN = fetNchans*nFetPerChan; % now number of features total
% fet = reshape(fet, size(fet,1), []);

N = numel(clu);
assert(fetNchans <= size(fet, 3) && size(fet, 1) == N , 'bad input(s)')

clusterIDs = unique(clu);
unitQuality = zeros(size(clusterIDs));
contaminationRate = zeros(size(clusterIDs));
noise_rate = zeros(size(clusterIDs));

fprintf('%12s\tQuality\tContam\tPrcnt1.5msViol\tGoodness\tNspikes\n', 'ID'); % comment to suppress printing out the intermediate results

spike_clusters = readNPY('spike_clusters.npy');
spike_times = double(readNPY('spike_times.npy'));

for c = 1:numel(clusterIDs)
    
    theseSp = clu==clusterIDs(c);
    n = sum(theseSp); % #spikes in this cluster
    if n < fetN || n >= N/2
        % cannot compute mahalanobis distance if less data points than
        % dimensions or if > 50% of all spikes are in this cluster
        unitQuality(c) = 0;
        contaminationRate(c) = NaN;
        continue
    end
    
    fetThisCluster = reshape(fet(theseSp,:,1:fetNchans), n, []);
    
    % now we need to find other spikes that exist on the same channels
    theseChans = fetInds(c,1:fetNchans);
%     for f = 1:fetNchans
%         thisChanInds = fetInds==theseChans(f);
%         [chanInds,clustWithThisChan] = find(thisChanInds');
%         chanInds = chanInds(clustWithThisChan~=c);
%         clustWithThisChan = clustWithThisChan(clustWithThisChan~=c);                
%         
%         otherSpikes = ismember(clu, clusterIDs(clustWithThisChan));
%         nOtherSpikes = sum(otherSpikes);
%         
%         fetOtherClusters = zeros(nOtherSpikes, nFetPerChan, fetNchans);
%         nInd = 1;              
%         
%         for t = 1:numel(clustWithThisChan)            
%             thisCfetInd = chanInds(t);
%             theseOtherSpikes = clu==clusterIDs(clustWithThisChan(t));
%             fetOtherClusters(nInd:nInd+sum(theseOtherSpikes)-1,:,f) = ...
%                 fet(theseOtherSpikes,:,thisCfetInd);
%             nInd = nInd+sum(theseOtherSpikes);
%         end
%         
%     end

    % for each other cluster, determine whether it has at least one of
    % those channels. If so, add its spikes, with its features put into the
    % correct places
    nInd = 1; fetOtherClusters = zeros(0,size(fet,2),fetNchans);
    for c2 = 1:numel(clusterIDs)
        if c2~=c
            chansC2Has = fetInds(c2,:);
            for f = 1:length(theseChans)
                
                if ismember(theseChans(f), chansC2Has)
                    
                    theseOtherSpikes = clu==clusterIDs(c2);
                    thisCfetInd = find(chansC2Has==theseChans(f),1);
                    fetOtherClusters(nInd:nInd+sum(theseOtherSpikes)-1,:,f) = ...
                        fet(theseOtherSpikes,:,thisCfetInd);
                end                
                
            end
            if any(ismember(chansC2Has, theseChans))
                nInd = nInd+sum(theseOtherSpikes);
            end
        end
    end
                   
    fetOtherClusters = reshape(fetOtherClusters, size(fetOtherClusters,1), []);
    
    [uQ, cR] = maskedClusterQualityCore(fetThisCluster, fetOtherClusters);
    
    unitQuality(c) = uQ;
    contaminationRate(c) = cR;
    
    % percent within 1.5 ms spikes (hard coded sampling rate of 30kHz)
    st=spike_times(spike_clusters==clusterIDs(c)-1);
    sp_isi=find(diff(st/30000)<.0015);
    sp = loadKSdir(pwd);
    
    index = find(sp.cids==clusterIDs(c)-1);
    if isempty(index)
        gnum = 0;
    else
    gnum = sp.cgs(index);
    end
    if gnum ==1 
        goodness= 'MUA';
    elseif gnum==2
        goodness = 'good';
    elseif gnum == 0
        goodness = 'noise';
    elseif gnum ==3 
        goodness = 'not sorted';
    else
        goodness = 'missing';
    end
    if isempty(sp_isi)
        num_sp_isi=0;
    else
        num_sp_isi=length(sp_isi);
    end
    total_sp=length(st);
    noise_rate(c)=(num_sp_isi/total_sp)*100;
    
    fprintf('cluster %3d: \t%6.1f\t%6.2f\t\t%6.2f\t\t%s\t\t%d\n', clusterIDs(c), unitQuality(c), contaminationRate(c),noise_rate(c), goodness,numel(st)); % comment to suppress printing out the intermediate results
    
    if uQ>1000
        keyboard;
    end
    
end
