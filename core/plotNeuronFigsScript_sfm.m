clear

SNR = []; 
SNR1 = [];
ksRoot = pwd;                                                               % Set dir to analyze

% loadPars.loadPCs = false;
% loadPars.excludeNoise = true;
sp = loadKSdir(ksRoot);                                                     % see loadKSdir for varargin

inclCID = sp.cids(sp.cgs == 2);
st = sp.st;
clu = sp.clu;

pcFeat = sp.pcFeat;
pcFeatInd = sp.pcFeatInd;
spikeAmps = sp.tempScalingAmps;
figDir = fullfile(ksRoot, 'figs'); 
if ~exist(figDir, 'dir') 
    mkdir(figDir); 
end

params.dataType = sp.dtype;
params.filename = sp.dat_path;
d = dir(fullfile(ksRoot, params.filename)); 


nSamp = d.bytes / 2 / sp.n_channels_dat;
params.dataSize = [sp.n_channels_dat nSamp];
params.chanMap = readNPY(fullfile(ksRoot, 'channel_map.npy'));
params.Fs = sp.sample_rate;
params.xcoords = sp.xcoords; 
params.ycoords = sp.ycoords;
params.plotDistance = 100;
params.nWFsToLoad = 1000;
params.nWFsToPlot = 100;
params.gain = 0.6/512/500 * 1e6;                                            % raw file units to uV
params.nPCsToPlot = 50000;
% params.highlightRegions = [132 1104];

%% extract median WFs (just once) %% otherwise load them
cd(ksRoot)
inclSP = ismember(clu, sp.cids(sp.cgs == 2));

if ~exist('medWFs.mat', 'file')
    medWFs = extractMedianWFs(clu(inclSP), st(inclSP), params.Fs, params.filename, params.dataType, params.dataSize, params.chanMap, params.gain);
    save(fullfile(ksRoot, 'medWFs.mat'), 'medWFs');
else
    load(fullfile(ksRoot, 'medWFs.mat'))
end

%% compute cluster quality stats (just once)

[cids, cgs, uQ, cR, isiV, noise_rate] = sqKilosort.computeAllMeasures(pwd);

save(fullfile(ksRoot, 'clusterQualityMetrics.mat'), 'cids','cgs', 'uQ', 'cR', 'isiV', 'noise_rate');

%%

sparsePCfeat = sparsePCs(pcFeat, pcFeatInd, sp.spikeTemplates, 2, 10);

%%
figDir = fullfile(sprintf('%s',pwd),'\figs');

theseISI = isiV(cgs == 2);
theseCR = cR(cgs == 2);
theseID = uQ(cgs == 2);
maxChans = [];
for q = 1:length(inclCID)
    clusterID = inclCID(q);
    stats.medWF = squeeze(medWFs(inclCID == clusterID,:,:))';
    stats.isiContamination = theseISI(inclCID == clusterID);
    stats.isoDistance = theseID(inclCID == clusterID);
    stats.mahalContamination = theseCR(inclCID == clusterID);
    [figHand, maxChan, snr, snr1] = neuronFig(clusterID, st, clu, sparsePCfeat, spikeAmps, stats, params);
    set(figHand, 'Position', [-1890         -59        1810        1031]);
    savefig(figHand, fullfile(figDir, sprintf('/cluster%d', clusterID)))
    close(figHand); 
    clear figHand
    maxChans(q) = maxChan(1);
    SNR(q) = snr;
    SNR1(q) = snr1;
end
save('maxChans.mat', 'maxChans') %save max spike channels
save('SNR.mat', 'SNR', 'SNR1') %save signal to noise calculations
plotAllMeasures(cgs, uQ, cR,isiV) % plot all clusters and their measures of quality

[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(pwd); %plot spikes across time to look at drift
figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths);

%plot basic frequency characterization
lfpD = pwd; % LFP file from spikeGLX specifically
fn1 = dir('*.dat');
lfpFilename = fn1.name;

lfpFs = 30e3;  % neuropixels phase3a
nChansInFile = 64;  % neuropixels phase3a, from spikeGLX

[lfpByChannel, allPowerEst, F, allPowerVar] = ...
    lfpBandPower(lfpFilename, lfpFs, nChansInFile, []);

chanMap = readNPY(fullfile(pwd, 'channel_map.npy'));
nC = length(chanMap);

allPowerEst = allPowerEst(:,chanMap+1)'; % now nChans x nFreq

% plot LFP power
dispRange = [0 100]; % Hz
marginalChans = [1:2:nC];
freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};

plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands);