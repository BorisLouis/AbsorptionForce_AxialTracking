clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = ['F:\LAIA\2024\SPT\00. 2D CAL - BORIS'];

%file info
file.path  = 'F:\LAIA\2024\SPT\20240808_noFiltering\pH 2';
file.ext   = '.ome.tif';
path2Cal = 'F:\LAIA\2024\SPT\00. 2D CAL - BORIS';
dimension = '3D';

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 80;
detectParam.consThresh = 6;
detectParam.FWHM_pix = 8;
%tracking parameter
trackParam.radius  = 2000;%nm
trackParam.memory  = 3;
pxSize = 95;
expTime = 0.01;% in second

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'MaxLR'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = 'true'; %true to recalibrate;

%% Calibrate data
mainDir = dir(file.path);
mainDir = mainDir(cell2mat({mainDir.isdir}));
%loop through the content of the directory
count = 0;
allTraces = {};
for i = 3:size(mainDir,1)
    %Check if the directory
    folderPath = [mainDir(i).folder filesep mainDir(i).name];
    file2Analyze = Core.Movie.getFileInPath(folderPath,file.ext);
    % Calibrate data
    if ~isempty(file2Analyze)
        count = count+1;
        f.path = file2Analyze.folder;
        f.ext = file.ext;
        currFile = Core.MPMovie(f , path2Cal,info);

        if count == 1
            currFile.giveInfo;
            allInfo = currFile.getInfo;
        else
            currFile.info = allInfo;
        end
        currFile.calibrate;
        
        
    else

        warning([mainDir(i).folder filesep mainDir(i).name ' did not contain any ' obj.ext ' file and is therefore ignored']);

    end
    
    % Localization
    nF = currFile.raw.maxFrame(1);
    
    position = table(zeros(500,1),zeros(500,1),zeros(500,1),zeros(500,1),...
                   'VariableNames',{'row', 'col','plane','frame'});
    for j = 1:8
        mov = double(currFile.getPlane(j));

        %getFrame with highest Value 
        int = mov(:);
        T = otsuthresh(int(:)./max(int(:)));
        bwImage = (mov./max(int(:)))<T/3;
        bwImage = imcomplement(bwImage);
        
        for k=1:size(bwImage,3)
            
            bw = bwareaopen(bwImage(:,:,k),300);
            SE = strel('disk',5);
            bw = imopen(bw,SE);
           
%             if k ==62
%                 disp('stop')
%             end
%             
            [ctr] = regionprops(bw,'Area','Centroid');          
            
            pos = cat(1,ctr.Centroid);
            if ~isempty(pos)
            startIdx = find(position.row==0,1,'First');
            if isempty(startIdx)
                startIdx = length(position.row)+1;
            end
            
            %add focus metric based on gradient, determine xy by gradFit
            %too
             pos = flip(pos,2);
             pos(:,3) = j;
             pos(:,4) = k;
             position(startIdx:startIdx+size(pos,1)-1,:) = array2table(pos);
            end
        end
    end
%% Consolidate position
    consPos = table(zeros(500,1),zeros(500,1),zeros(500,1),zeros(500,1),zeros(500,1),...
                   'VariableNames',{'row', 'col','frame','plane','ID'});
    for j = 1:nF
       
        pos = [position.row(position.frame==j), position.col(position.frame==j), position.plane(position.frame==j)];
%         if j == 62
%            disp('stop') 
%         end

        
        ct = 0;
        while ~isempty(pos)
            ct = ct +1;
            currPos = pos(1,:);
            
            match = pos(and(abs(pos(:,1)-currPos(1))<8,abs(pos(:,2)-currPos(2))<8),:);
            
            match(:,4) = match(:,3);
            match(:,3) = j;
            match(:,5) = ct;
            if length(unique(match(:,4)))>3
                startIdx = find(consPos.row==0,1,'First');
                
                if isempty(startIdx)
                    startIdx =height(consPos);
                end
                    consPos(startIdx:startIdx+size(match,1)-1,:) = array2table(match);
    
            else
                
            end
            %remove from list
            idx = ismember(pos,match);
            pos(idx(:,1),:) = [];
            
        end
    end
    consPos(consPos.row==0,:) = [];
%% 3D info
    %get 3D info
    nJ = height(consPos);
    list = consPos;
    SRList = table(zeros(500,1),zeros(500,1),zeros(500,1),zeros(500,1),...
                   'VariableNames',{'row', 'col','z','t'});
    ct = 0;
    while ~isempty(list)
        ct = ct+1;
        currPart = list(and(list.ID == list.ID(1),list.frame == list.frame(1)),:); 
        
        volIm = currFile.getFrame(currPart.frame(1));
        
%         if list.frame(1)==58
%             disp('stop');
%             
%         end
%         
%         if list.frame(1)==87
%             disp('stop');
%         end

        SRList.row(ct) = mean(currPart.row);
        SRList.col(ct) = mean(currPart.col);
        SRList.t(ct) = currPart.frame(1);
        
        [partVolIm] = getPartVolIm(SRList(ct,:),15,volIm);
        
        if SRList.t(ct) == 165
            disp('stop')
        end
        
        if ~isempty(partVolIm)
            planePos = currFile.calibrated.oRelZPos;
            SRList.z(ct) = calcZ(partVolIm,planePos);
            
%             if ~isnan(SRList.z(ct))
%                 [val,idx] = min(abs(SRList.z(ct)-planePos));
%                 [x,y,~,~] = Localization.gradFit(partVolIm(:,:,idx),round(size(partVolIm,1)/2)-3);
%             else
%                 x = NaN;
%                 y = NaN;
%             end
 
            %previously
            SRList.row(ct) = mean(currPart.row)*pxSize;
            SRList.col(ct) = mean(currPart.col)*pxSize;
%             SRList.row(ct) = round(size(partVolIm,1)/2)+y*pxSize;
%             SRList.col(ct) = round(size(partVolIm,1)/2)+x*pxSize;
            SRList.t(ct) = currPart.frame(1);
            %disp('stop')
        else
            
            SRlist.row(ct) = 0;
            SRList.col(ct) = 0;
            SRList.t(ct) = currPart.frame(1);
            
        end
        
        %get Image
        %get part Volume
        %get 3D position
        list(and(list.ID == list.ID(1),list.frame == list.frame(1)),:) = [];
        
        
    end
    % clean the list
   
    %%
    SRList(SRList.row==0,:) = [];
    SRList(isnan(SRList.z),:) = [];
    %tracking
    traces = trackParticle(currFile,SRList,trackParam);
   
    traces(cellfun(@height,traces)<=10) = [];
    
    fileID = ones(length(traces),1)*i-2;
    fileID = num2cell(fileID);
    traces = traces';
  
    %% calculate speed
    speedXY = zeros(length(traces),1);
    speedZ = zeros(length(traces),1);
    deltaZ = zeros(length(traces),1);
    
    

    for j = 1:length(traces)
       currTrace = traces{j};
       deltaZ(j) = max(currTrace.z)-min(currTrace.z);
       %we only calculate for particles that were followed along
       %significant z
       
       speedXY(j) = mean([mean(diff(currTrace.row)),mean(diff(currTrace.col))]);
       speedZ(j)  = mean(diff(currTrace.z));
           
       

    end  
    
      allTraces = [allTraces; cat(2, traces, fileID, num2cell(deltaZ(:)), num2cell(speedXY),num2cell(speedZ))];
    
end


%% Plot
figure
hold on
for i =1:size(allTraces,1)
    
    currTrace = allTraces{i,1};
  
    colPlot = currTrace.col;
    rowPlot = currTrace.row;
    zPlot   = currTrace.z;
    tPlot   = currTrace.rT-currTrace.rT(1);
    plot3(colPlot,rowPlot,zPlot)
    %plot with time color coding
    patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')
   
    
end
zlim([-2200 2200])
colorbar
caxis([0 4000])
box on

%% calculate mean speed and std and unit conversion
%nm/exposure to micrometer per second
%structure of the data:
%ALLDATA(:,1)
%contains the tracked data for a specifc trace, containing the
%x,y,z,frame,real time and particle ID
%ALLDATA(:,2)
%ID of the movie from which the traces was detected
%AllData(:,3)
%Speed in XY (average X, average Y and then average of both)
%AllData(:,4)
%Speed in Z (average for the trace)

meanSXY = mean(cell2mat(allTraces(:,4)))/1000/expTime;
stdXY   = std(cell2mat(allTraces(:,4)))/1000/expTime;
meanSZ = mean(cell2mat(allTraces(:,5)))/1000/expTime;
stdZ   = std(cell2mat(allTraces(:,5)))/1000/expTime;

figure
hold on
errorbar(meanSXY,stdXY,'o')
errorbar(meanSZ,stdZ,'o')
legend({'XY','Z'})
box on
axis square

%% Mean squared Displacement
msdXY = zeros(size(allTraces,1),max(SRList.t));
msdZ  = msdXY;
for i = 1:size(allTraces,1)
    currPart = allTraces{i,1};
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


plottingRange = find(isnan(mean(msdXY,1)),1,'first')-1;

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

MSD_all = cat(1,MSDXY(1:plottingRange),MSDZ(1:plottingRange),sMSDXY(1:plottingRange),sMSDZ(1:plottingRange));
filename = [file.path filesep 'MSD.mat'];
save(filename,'MSD_all')

disp('========= Analysis summary ========')
disp(['Total number of particles tracked: ' num2str(size(allTraces,1))])
disp(['Average speed in XY :' num2str(round(meanSXY)) ' um/sec'])
disp(['Standard deviation speed in XY: ' num2str(round(stdXY))])
disp(['Average Speed in Z: ' num2str(round(meanSZ)) ' um/sec'])
disp(['Standard deviation speed in Z: ' num2str(round(stdZ))])

disp('======= End Analysis summary =======')



function [partVolIm] = getPartVolIm(partData,ROIRad,volIm)
            %Extract data from particle in the 8 planes
    imSize = size(volIm);

    pos = [round(nanmean(partData.col)),round(nanmean(partData.row))];

    ROIs = Misc.getROIs(pos,ROIRad,imSize(1:2));
    
    try
        
        partVolIm = volIm(ROIs(3):ROIs(4),ROIs(1):ROIs(2),:);
    catch e
        partVolIm = [];
    end

end

function z = calcZ(partVolIm,planePos)

for i = 1:size(partVolIm,3)
   fMet(i,:) =  mean(partVolIm(:,:,i),1) + mean(partVolIm(:,:,i),2)';
   testMetc(i,:) = max(mean(imgradient(partVolIm(:,:,i))));
end

%Get ROI XZ, YZ scaled to same pixel size
% [Mag] = Core.MPLocMovie.getZPhasorMag(partVolIm);

domain = planePos;
%data   = [Mag.x]+[Mag.y];
data   =  testMetc';
guess.sig = 1;
guess.mu  = domain(testMetc==max(testMetc));
guess.A   = max(data)-min(data);

%             [Res,fitData] = SimpleFitting.gauss1D(data,domain,guess);
%             RMS = sqrt(sum((data(:)-fitData(:)).^2)/length(data));
%             adjR = 1 - RMS.^2/var(data);
%           
params = [guess.sig guess.mu guess.A min(data)];

fun = @(x) SimpleFitting.minGauss1D(domain,data,x);
opt = optimset('Display','off');
% then we can do:
[out, RMSD] = fminsearch(fun,params,opt);
%normalize RMSD with mean data
adjR = 1 - RMSD.^2/var(data);

%[~, gaussFit] = SimpleFitting.minGauss1D(domain,data,out);

z = out(2);
%if the z position is out of bound we do not consider the data
if or(z<min(domain)*1.2,z>max(domain)*1.2)
    z   = NaN;                           
  
else
    z = z*1000;

end

end

function traces = trackParticle(file, SRList,trackParam)
          
            DataToTrack = SRList;
            ImMax = max(DataToTrack.t);
            % get the timing of each frame
            if ~isfield(file.raw.movInfo,'timing')
                timing = ((1:ImMax)-1)*file.raw.movInfo.expT;
            else
                timing = file.raw.movInfo.timing;
            end
            allTimes = timing(DataToTrack.t);
            DataToTrack.rT = allTimes(:);
            %Converts data
            [ToTrack,AllField] = Core.trackingMethod.ConvertData(DataToTrack,ImMax);
            count = 0; 
            %check that there are particles in frame 1
            if or(isempty (ToTrack{1}),isempty(ToTrack{2}))
                ToTrack(1) = [];
                count = 1;
                while or(isempty(ToTrack{1}),isempty(ToTrack{2}))
                    ToTrack(1) = [];
                    count = count+1;
                end   
            end    
            
            Initialized = [ToTrack{1},(1:size(ToTrack{1},1))'];

            ToTrack(1) = [];

            %%%%% INITIALIZE FOR TRACKING DATA
            % TrackedData = Tracked;
            if count >0
                TrackedData_data = cell(1,count);
                TrackedData_data = [TrackedData_data,{Initialized}];
            else
                TrackedData_data = {Initialized};
            end
            
            TrackedData_maxid = max(TrackedData_data{end}(:,end));

            NextFrame = Core.Tracking.CoordsInFrameNextFrame;
            NextFrame.dataNext = ToTrack{1};
            NextFrame.timeNext = NextFrame.dataNext(1,end);

            MemoryArray_data = [];

            %%%%% TRACK DATA RECURSIVELY
            radius = trackParam.radius;
            MaximumTimeMem = trackParam.memory;
            totlengthFinal =0;
            h = waitbar(0,'Tracking particles...');
            while ~isempty(ToTrack) 
                if ~isempty(NextFrame.dataNext)
                    waitbar(NextFrame.timeNext/ImMax,h,'TrackingParticles');
                    %Upon every iteration, memory must be cleaned to remove any particles
                    %that have possible moved out of the frame, and have stayed out of the
                    %frame for longer then MaximumTimeMem

                    MemoryArray_data = Core.trackingMethod.CleanMemory(MemoryArray_data,NextFrame.timeNext,MaximumTimeMem);

                    %Add remaining memory to the previous frame and initialize the array
                    %for tracking

                    ImPrev = Core.Tracking.CoordsInFramePreviousFrame;

                    ImPrev.time = length(TrackedData_data);
                    ImPrev.data = ImPrev.AddMemoryToPreviouslyTrackedData( TrackedData_data{end},MemoryArray_data);

                    %Search for neighbours in the immedeate vicitinity of radius = radius
                    NextFrame.PossibleNeighbours = Core.trackingMethod.SearchNeighbours(NextFrame,ImPrev.data,radius);

                    %Resolve any conflitcs by minimizing the overall distance. 

                    NextFrame.PossibleNeighbours = Core.trackingMethod.ResolveConflicts(NextFrame.PossibleNeighbours);

                    %Add any unassigned particles to the memory

                    MemoryArray_data = Core.trackingMethod.AddToMemory(NextFrame.PossibleNeighbours, ImPrev.data);

                    %Add tracked data to the array of tracked data by assigning possible
                    %candidates to the indeces of the next frame

                    [TrackedData_data{end+1},TrackedData_maxid] = Core.trackingMethod.AddDataToTracked( ImPrev.data ,NextFrame.dataNext, NextFrame.PossibleNeighbours,TrackedData_maxid);

                end
                ToTrack(1) = [];
                
                %Boris if the tracked data is empty we store the data f
                if and(isempty(TrackedData_data{end}),~isempty(NextFrame.dataNext))
                    Initialized = [NextFrame.dataNext,(TrackedData_maxid+1:TrackedData_maxid+size(NextFrame.dataNext,1))'];
                    TrackedData_data{end} = Initialized;
                    TrackedData_maxid = max(TrackedData_data{end}(:,end));
                end
                
                %Initialize for next step
                if ~isempty(ToTrack)
                NextFrame.dataNext = ToTrack{1};
                    if ~isempty(ToTrack{1})
                        NextFrame.timeNext =NextFrame.dataNext(1,end);
                      
                    else
                        NextFrame.timeNext =NextFrame.timeNext+1;
                        %Boris Fix for data with "holes"
                        TrackedData_data{end+1}= [];
                        
                    end
                end

            end
            close(h);
            AllParticles = cell(1,TrackedData_maxid);
            
            TrackedData = Core.trackingMethod.ConvertFinalOutput( TrackedData_data,AllParticles,AllField);
            
%             file.particles.traces = TrackedData;
%             file.particles.nTraces = length(TrackedData);
            
            %[trace3D] = obj.get3DTraces;
            traces = TrackedData;
%             file.traces3D = TrackedData;
%             
%             filename =[file.raw.movInfo.Path filesep 'Traces3D.mat'];
%             
%             save(filename,'TrackedData');
%             
            
end
        


