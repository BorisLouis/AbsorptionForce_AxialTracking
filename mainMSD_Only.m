
expTime = mean(diff(traces2Keep{1}.rT))/1000; %extract exposure time

meanSXY = mean(cell2mat(traces2Keep(:,4)))/1000/expTime;
stdXY   = std(cell2mat(traces2Keep(:,4)))/1000/expTime;
meanSZ = mean(cell2mat(traces2Keep(:,5)))/1000/expTime;
stdZ   = std(cell2mat(traces2Keep(:,5)))/1000/expTime;

figure
hold on
errorbar(meanSXY,stdXY,'o')
errorbar(meanSZ,stdZ,'o')
legend({'XY','Z'})
box on
axis square

%% Mean squared Displacement
msdXY = zeros(size(traces2Keep,1),max(SRList.t));
msdZ  = msdXY;
for i = 1:size(traces2Keep,1)
    currPart = traces2Keep{i,1};
    coordXY = [currPart.row  currPart.col];
    coordZ  = currPart.z;
    
    msdXY(i,1:length(coordXY)-1) = MSD.calc(coordXY/10^3);
    msdZ(i,1:length(coordZ)-1) = MSD.calc(coordZ/10^3);
    
end

msdXY(msdXY==0) = NaN;
msdZ(msdZ ==0)  = NaN;
MSDXY = nanmean(msdXY,1);
MSDZ  = nanmean(msdZ,1);
sMSDXY = nanstd(msdXY,1);
sMSDZ  = nanstd(msdZ,1);

xAxMSD = (1:length(MSDXY))*expTime;


plottingRange = 3*find(isnan(mean(msdXY,1)),1,'first')-1;

figure
hold on
h1 = Plotting.shadedErrorBar(xAxMSD(1:plottingRange),MSDXY(1:plottingRange)/2,sMSDXY(1:plottingRange)/2,'lineProps','-b');
h2 = Plotting.shadedErrorBar(xAxMSD(1:plottingRange),MSDZ(1:plottingRange),sMSDZ(1:plottingRange),'lineProps','-r');
xlabel('Lag time (s)')
ylabel('MSD(\mum^2/s)')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
axis square
box on
legend([h1.mainLine h2.mainLine],'XY','Z')

MSD_all = cat(1,MSDXY(~isnan(MSDXY)),MSDZ(~isnan(MSDXY)),sMSDXY(~isnan(MSDXY)),sMSDZ(~isnan(MSDXY)));
filename = [file.path filesep 'MSD.mat'];
save(filename,'MSD_all')