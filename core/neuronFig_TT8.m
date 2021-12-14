function [Hfig, maxChan] = neuronFig_TT8(clusterID, spikeTimes, clu, sparsePCfeat, spikeAmps, stats, params)
% Makes a plot of relevant stats/figures for a neuron
%
% stats has:
% clusterID
% medWF (uVolts)
% snr,snr1
% - isoDistance
% - isiContamination
% - mahalContamination
%
% params has:
% - chanMap
% - xcoords, ycoords
% - "dataType" and "dataSize" of raw file
% - raw "filename"
% - Fs

% 
% if isempty(stats)
%     stats.isiContamination = NaN;
%     stats.isoDistance = NaN;
%     stats.snr1 = nan;
%     stats.snr = nan;
%     stats.medWF = nan;
% end

Fs = params.Fs;
wfWin = -round(1*Fs/1000):round(2.5*Fs/1000); nWFsamps = numel(wfWin);
nChInFile = params.dataSize(1);
nSamp = params.dataSize(2);
nCh = numel(params.chanMap);
xc = params.xcoords; yc = params.ycoords;

theseST = spikeTimes(clu==clusterID);
%params.nWFsToLoad = params.nWFsToPlot;
nWFsToLoad = min(params.nWFsToLoad, length(theseST));
extractST = round(theseST(randperm(length(theseST), nWFsToLoad))*Fs);
mmf = memmapfile(params.filename, 'Format', {params.dataType, [nChInFile nSamp], 'x'});

theseWF = zeros(nWFsToLoad, nCh, nWFsamps);
for i=1:nWFsToLoad
    tempWF = mmf.Data.x(1:nChInFile,extractST(i)+wfWin(1):extractST(i)+wfWin(end));
    theseWF(i,:,:) = tempWF(params.chanMap+1,:) * params.gain;
end

%medWF = squeeze(median(double(theseWF),1));
%medWFuV = medWF.*gain;
medWFuV = stats.medWF;

chanRMS = sqrt(mean(medWFuV.^2,1));
maxChan = find(chanRMS==max(chanRMS),1);

Hfig = figure;
set(Hfig, 'Color', 'w');
set(Hfig,'name',['cluster# ' num2str(clusterID)])
neuronColor = [ 0    0.4470    0.7410];
otherColor = [0.4660    0.6740    0.1880];

% plot of waveforms
Hax = axes('position', [.05 .35 .3 .6]); hold on;
plot((double(squeeze(theseWF(:,maxChan,:))) - mean(squeeze(theseWF(:,maxChan,:)),2))','color',neuronColor)
medWF = squeeze(median(double(theseWF),1));
plot(medWF(maxChan,:)-mean(medWF(maxChan,:)),'color',otherColor,'linewidth',2)
ylim([min(medWF(maxChan,:)-mean(medWF(maxChan,:))) max(medWF(maxChan,:)-mean(medWF(maxChan,:)))]*6)
xlim([0 nWFsamps]);
title('waveforms')
box off;

% compute PC stuff
thesePCs = sparsePCfeat(clu==clusterID,:);
meanPC = mean(thesePCs);
[~, ii] = sort(abs(meanPC), 2, 'descend');
topChans = ii(1:2);

% Method 1: just pick the top two channels for this cluster
otherSpikesIncl = sparsePCfeat(:,topChans(1))~=0 & ...
    sparsePCfeat(:,topChans(2))~=0 & clu~=clusterID;
% otherSpikesPCs = sparsePCfeat(otherSpikesIncl, topChans);
% np = min(size(otherSpikesPCs,1), params.nPCsToPlot);
% otherPCsToPlotInds = randperm(size(otherSpikesPCs,1), np);
% otherPCsToPlot = otherSpikesPCs(otherPCsToPlotInds,:);
thesePCsToPlot = thesePCs(:,topChans);

% plot of PC space
Hax = axes('position', [.40 .55 .25 .3]);
hold on;
plot(thesePCsToPlot(:,1), thesePCsToPlot(:,2),'.', 'MarkerSize', 3, 'Color', neuronColor);
drawnow;
title({['PCs iso dist = ' num2str(round(stats.isoDistance,2))];['   Mahal contam = ' num2str(round(stats.mahalContamination,2))]})
box off

% plot of stability over time
Hax = axes('position', [.05 .1 .45 .15]); hold off;
if isfield(params, 'highlightRegions') && ~isempty(params.highlightRegions)
    % field is nx2
    for q = 1:size(params.highlightRegions,1)
        hr = params.highlightRegions(1,:);
        fill([hr(1) hr(2) hr(2) hr(1)], [0 0 1000 1000], [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
        hold on;
    end
end

binSize = 20;
timeBins = 0:binSize:ceil(spikeTimes(end));
[n,x] = hist(theseST, timeBins);
n = n./binSize;

tpc = thesePCs(:,topChans(1));
tpc = (tpc-mean(tpc))./std(tpc)*std(n);
tpc = tpc-min(tpc);

h1 = plot(theseST, tpc, '.', 'Color', neuronColor);
hold on;

h2 = stairs(x,n, 'LineWidth', 2.0, 'Color', 'k');

ymax = max([max(n) max(tpc)]);
ylim([0 ymax]);

title('top pc, firing rate over time')
xlabel('time(sec)')
H = ylabel('firing rate (sp/sec)');
pos = get(H,'position');
set(H,'position',[pos(1)*.9 pos(2) pos(3)])
H = legend([h1 h2], {'top PC', 'f rate'});
set(H,'position',[.51 .13 .10 .07])
box off;

% histogram of fit amplitudes
Hax = axes('position', [.78 .1 .15 .18]);

h = histogram(spikeAmps(clu==clusterID), 'Normalization', 'pdf', ...
    'DisplayStyle', 'stairs');
n = get(h, 'Values'); x = get(h, 'BinEdges');
stairs(x(1:end-1), n, 'LineWidth', 2.0);

xlabel('relative spike amplitude')
ylabel('pdf')
title('PDF spike amplitude');

% histogram of ISIs
Hax = axes('position', [.72 .55 .25 .3]);
ISI = diff(theseST)*1000; % milliseconds
ISIviolations2 = (length(find(ISI<2))/length(ISI))*100;
ISIviolations15 = (length(find(ISI<1.5))/length(ISI))*100;
bins = 0.001:.1:2000;
H = histc(ISI,bins);
semilogx(bins,H)
YL = ylim;
hold on
semilogx([2 2],YL,'k:','linewidth',1)
xlim([1 1000])
xlabel('log10(ISI) (ms)')
ylabel('# spikes')
title({'ISI histogram'; ['<2ms: ' num2str(round(ISIviolations2,2)) '%']; ['<1.5ms: ' num2str(round(ISIviolations15,2)) '%']});
box off

% misc stats printed here
Hax = axes('position', [.40 .32 .30 .1]);
text(0,1,{['nSpikes: ' num2str(length(theseST))]; ['SNR: ' num2str(round(stats.snr1,2))];['SNR(rqq): ' num2str(round(stats.snr,2))]})
axis off

