% plotKSdata_TT8
% in progress 5/16/19 for a quick listing of quality etc. for 'good' clusters

flag.force_calc_medWFs = 1;
flag.force_calc_AllMeasures = 1;
flag.draw_clusterFig = 1;

ksRoot=pwd;

Fs = 30000;
nChInFile = 32;

loadPars.loadPCs = true;
sp = loadKSdir(ksRoot, loadPars);

inclCID = sp.cids(sp.cgs==2);
st = sp.st;
clu = sp.clu;

pcFeat = sp.pcFeat;
pcFeatInd = sp.pcFeatInd;
spikeAmps = sp.tempScalingAmps;
figDir = fullfile(ksRoot, 'figs');
if ~exist(figDir,'dir'); mkdir(figDir); end

params.dataType = sp.dtype;
params.filename = sp.dat_path;
d = dir(fullfile(ksRoot,params.filename));

nSamp = d.bytes/2/sp.n_channels_dat;
params.dataSize = [sp.n_channels_dat nSamp];
params.chanMap = readNPY(fullfile(ksRoot, 'channel_map.npy'));
params.Fs = sp.sample_rate;
params.xcoords = sp.xcoords; params.ycoords = sp.ycoords;
params.plotDistance = 100;
params.nWFsToLoad = 1000;
params.nWFsToPlot = 100;
params.gain = 0.6/512/500*1e6; % raw file units to uV
params.nPCsToPlot = 50000;

% compute (or load) cluster quality stats
if ~exist(fullfile(ksRoot, 'clusterQualityMetrics.mat'),'file')
    flag.force_calc_AllMeasures =1;
end
if flag.force_calc_AllMeasures
    clusterPath = fullfile(pwd, 'cluster_groups.csv');
    spikeClustersPath = fullfile(pwd,'spike_clusters.npy');
    spikeTemplatesPath = fullfile(pwd,'spike_templates.npy');
    
    if exist(clusterPath, 'file')
        [~, cgs] = readClusterGroupsCSV(clusterPath);
    elseif exist(spikeClustersPath, 'file')
        clu = readNPY(spikeClustersPath);
        cgs = 3*ones(size(unique(clu)));        % all unsorted
    else
        clu = readNPY(spikeTemplatesPath);
        cgs = 3*ones(size(unique(clu)));        % all unsorted
    end
    % cgs: 0) noise, 1) , 2) good
    [cids, uQ, cR, noise_rate] = sqKilosort.maskedClusterQuality(pwd);
    % cids: clusterID
    % uQ: quality
    % cR: Mahal contamination
    % noise_rate: percent 1.5 ms ISI violations    
    save(fullfile(ksRoot, 'clusterQualityMetrics.mat'), 'cids','cgs', 'uQ', 'cR', 'noise_rate');
else
    load(fullfile(ksRoot, 'clusterQualityMetrics.mat'))
end

% extract (or load) median WFs
if ~exist(fullfile(ksRoot, 'medWFs.mat'),'file')
    flag.force_calc_medWFs = 1;
end
if flag.force_calc_medWFs
    cd(ksRoot)
    inclSP = ismember(clu, sp.cids(sp.cgs==2));
    medWFs = extractMedianWFs(clu(inclSP), st(inclSP), params.Fs, params.filename, ...
        params.dataType, params.dataSize, params.chanMap, params.gain);
    save(fullfile(ksRoot, 'medWFs.mat'), 'medWFs');
else
    load(fullfile(ksRoot, 'medWFs.mat'))
end

wfWin = -round(1/1000*Fs):round(3.5/1000*Fs);
nWFsamps = numel(wfWin);
%theseISI = isiV(cgs==2);
theseCR = cR(cgs==2);
theseID = uQ(cgs==2);
maxChans=nan(1,length(inclCID));
SNR = zeros(1,length(inclCID));
SNR1 = zeros(1,length(inclCID));
    
for iclu = 1:length(inclCID)
    clusterID = inclCID(iclu);
    theseST = st(clu==clusterID);
    stats.clusterID = clusterID;
    
    % Extract waveforms from this neuron and from background
    nWFsToLoad = min(params.nWFsToLoad, length(theseST));
    extractST = round(theseST(randperm(length(theseST), nWFsToLoad))*Fs);
    
    extractBckg = randperm(nSamp-numel(wfWin)-2, params.nWFsToPlot)-wfWin(1); % *samples* of background spikes
    mmf = memmapfile(params.filename, 'Format', {params.dataType, [nChInFile nSamp], 'x'});
    
%    theseWF = zeros(nWFsToLoad, nChInFile, nWFsamps);
    bckgWF = zeros(params.nWFsToPlot, nChInFile, nWFsamps);
%     for i=1:nWFsToLoad
%         % takes time:
%         tempWF = mmf.Data.x(1:nChInFile,extractST(i)+wfWin(1):extractST(i)+wfWin(end));
%         theseWF(i,:,:) = tempWF(params.chanMap+1,:);
%     end
    
    for i = 1:params.nWFsToPlot
        % takes time:
        tempWF = mmf.Data.x(1:nChInFile,extractBckg(i)+wfWin(1):extractBckg(i)+wfWin(end));
        bckgWF(i,:,:) = tempWF(params.chanMap+1,:);
    end
    
    % need median waveforms here
    stats.medWF = squeeze(medWFs(inclCID==clusterID,:,:))';
    medWFuV = stats.medWF;
    
    % choose channels to plot
    chanAmps = max(medWFuV)-min(medWFuV);
    maxChan = find(chanAmps==max(chanAmps),1);
    
    wfAmp = max(chanAmps);
    
    bckgMaxChan = double(squeeze(bckgWF(:,maxChan,:)));
    stats.snr1 = wfAmp./std(bckgMaxChan(:));
    stats.snr = wfAmp./median(abs(bckgMaxChan(:))/0.6745); % RQQ method
    
    %stats.isiContamination = theseISI(inclCID==clusterID);
    stats.isoDistance = theseID(inclCID==clusterID);
    stats.mahalContamination = theseCR(inclCID==clusterID);
    if flag.draw_clusterFig
        sparsePCfeat = sparsePCs(pcFeat, pcFeatInd, sp.spikeTemplates, 2, 10);
        figDir = fullfile(sprintf('%s',pwd),'\figs');
        [figHand, maxChan] = neuronFig_TT8(clusterID, st, clu, sparsePCfeat, spikeAmps, stats, params);
        %set(figHand, 'Position', [-1890         -59        1810        1031]);
        savefig(figHand, fullfile(figDir, sprintf('/cluster%d', clusterID)))
    end
    %     close(figHand); clear figHand
    maxChans(iclu)=maxChan(1);
    SNR(iclu)=stats.snr;
    SNR1(iclu)=stats.snr1;
    
end

% save('maxChans.mat', 'maxChans') %save max spike channels
% save('SNR.mat', 'SNR', 'SNR1') %save signal to noise calculations
% plotAllMeasures(cgs, uQ, cR,isiV) % plot all clusters and their measures of quality


