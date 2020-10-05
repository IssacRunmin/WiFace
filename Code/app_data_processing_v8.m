function Re = app_data_processing_v8(app)
global Running Log FNum AxesHandle LogStClock;
% Use in WiBlinkAuth.mlapp
% CSI & magnitude change:
% ------------------------------+++
% CSI Read:                     |-|
% mag :                  1234...
% -----------------------+++++++++++++++++++++***
% | mag did not contain  | | change_plot_time | new data
%                        | plot_start_time
% changed mag:             1234...
% -------------------------++++++++++++++++++++++
%                          | plot_start_time
% 1. Intialize Variables/ parameters
% repeat
% 2. Read the file 'com' to get the file name
% 3. Loop the following codes until the file didn't be written
%   3.1 Copy the file and Read the file for the CSI raw data
%   3.2 Data processing
%   3.3 Draw figure
% go to 2.

% 1. Initialize variables
% Data processing 
    Process = {
        'MultipathRemoval', 0
        'CSI_Ratio', 0
        'OutlierRemoval', 1
        'PhaseCalculation', 0
        'PCA', 0
        'BandpassFilter', 0
        'PowerDensitySpectrum', 0
        'FrequencyDomain', 0
        'PowerSpectralDensity', 0
        'AccelerationEstimation', 0
        'MeanAbsoluteDeviation', 0
        'DiscreteWaveletTransform', 0
        'LogTotalTime', 1
        'CountAction', 1
    };
% Process parmeters
    % CSI Parameters
    NSub = 30;
    NAntenna = 3;
    % Multipath removal
    MRCutTime = 0.5;
    % Bandpass filter
    CutoffFreq = [2 6];%[0.3 3];
    FilterOrder = 2;
    CutTime = 0.4; % second(s)
    % Hampel Identifiler / Outlier Removal
    HampelAntenna = 2;
    K0 = 600;
    Nsigma0 = 1.30;%1.30
    K1 = 200; % 200
    Nsigma1 = 1.70;%1.30
    % PCA
    NumPCA = 3;
    % Acceleration Estimation
    DownNum = 8;
%     Tau = 4e5;
    MinDist = 0.3; % second(s)
    % DWT - Discrete Wavelet Transform
    DWTlevel = 5;
    DWTwavelet = 'sym4';
    % Mean Absolute Deviation
    MADwindow = 1000;
    % Segmentation
    ActDet = [3 8 3 9]; 
        % (Feq1, Feq2, Threshold, DetectStartTime)
    WalkTime = [2, 8]; % minimum and maximum time of walking event
        % the minimum time is the time window.
    
% For Showing Data
    DebugTime = false;      % debug
    TimeOver = 5;      % Time out threshold: s
    MaxPlot = 10;  % second(s), Figure's x axis
    RefreshInterval = 0.5; % second(s), Refresh interval
% communication file
    CommFile = '/home/issacrunmin/MATLAB/com';
    ModelFile = 'home/issacrunmin/MATLAB/Mdl.mat';
% Error if file has already started
%% Initialization
HaveClassifier = true;
try
    load(ModelFile, 'Classifier', 'ConfidenceThr', 'Users');
catch ME
    if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
%         disp(ME.message);
        warning('Orm: Do not have model!');
        HaveClassifier = false;
    end
end
Process(:,2) = cellfun(@(x) {x==1},Process(:,2));
opts = cell2struct(Process(:,2),Process(:,1),1);
CommFid = fopen(CommFile,'r');
DataFileName = fgetl(CommFid);
fclose(CommFid);
state = DataFileName(end);
if state == '1'
    error('ERROR: Please Run This Script First!');
end
disp('Ready to Begin!');
while Running
    while state == '1'
        CommFid = fopen(CommFile,'r');
        DataFileName = fgetl(CommFid);
        fclose(CommFid);
        state = DataFileName(end);
    end
    while state ~= '1' && Running
        CommFid = fopen(CommFile,'r');
        DataFileName = fgetl(CommFid);
        fclose(CommFid);
        state = DataFileName(end);
        pause(RefreshInterval);
    end
% intialize
    DataFileName(end - 1 : end) = [];
    DataFileName(1 : 7) = [];
%     Duration = RefreshInterval;        % refresh time
    PlotSt = 1;            % used in plot()
    TotalTime = 0;                 % total time
    CursorSeries = ones(1,6000);
    Mag = zeros(NSub, 0, NAntenna);
    LogRaw = zeros(NSub, 0, NAntenna);
    Phase1 = zeros(NSub, 0, NAntenna);
    Phase2 = zeros(NSub, 0, NAntenna);
    RatioData = zeros(NSub, 0, NAntenna);
    LenAll = 0;
    LogStClock = clock;
    % File cursor
    CursorCurrent = 1;
    % Time Out
    TimeOutDur = 0;
    isBlinking = false;
    Slience = true;
    BlinkStarted = false;
    LogSt = TotalTime;
    tic
    while(TimeOutDur < TimeOver) && Running
        LastTotalTime = TotalTime;
        TotalTime = toc; % TotalTime + Duration;
        CursorLast = CursorCurrent;
%         tic
    % 3.1 Copy the file and Read the file for the CSI raw data
        % copyfile(file_name,'temp.dat','f');
        % Do not copy the file
%         try 
            [CSI_trace,CursorCurrent] = read_bf_file_online(...
                DataFileName, CursorLast);
%         catch ME
%             warning('Error in read_bfee: Wrong beamforming matrix size.');
% %             rethrow(ME);
%         end
        if CursorLast == CursorCurrent
            Duration = toc - TotalTime;
            if (Duration <= RefreshInterval)
                pause(RefreshInterval-Duration); 
                Duration = RefreshInterval;
            end
            TimeOutDur = TimeOutDur + Duration;
            continue;
        else
            TimeOutDur = 0;
        end
        LenNew = size(CSI_trace, 1);
        LenAll = LenAll + LenNew;
        CursorSeries(ceil(TotalTime)) = LenAll;
        RawData = zeros(30,LenNew, 3);
        for i = 1:LenNew
            try
            CSIEntry = CSI_trace{i}; %for every packet
            csi = get_scaled_csi(CSIEntry); 
            % csi(Ntx, Nrx, Channel) 2 * 3 * 30
            csi = squeeze(csi(1, :, :))'; 
            RawData(:, i, :) = csi;
            catch ME
                warning('Error in get_scaled_csi()');
                disp(getReport(ME));
                continue;
            end
        end
%% Logging
        if Log
            if isempty(LogRaw)
%                 LogSt = TotalTime;
%                 LogStClock = clock;
%                 if opts.CountAction
%                     Groundtruth(app, -1);
%                 end
                Lag = etime(clock, LogStClock);
                LogSt = toc - Lag;
                StIdx = floor(Lag / (TotalTime - LastTotalTime) * LenNew);
                StIdx = min(max(StIdx, 1), LenNew);
                LogRaw = RawData(:,StIdx:LenNew,:);
            else
                LogRaw = cat(2, LogRaw,RawData);
            end
        else
            if ~isempty(LogRaw)
                if opts.CountAction
                    GroundTruth = Groundtruth(app,0);
                else
                    GroundTruth = [];
                end
                app.CountButton.Enable = 'off';
                app.DeleteButton.Enable = 'off';
                app.StartButton.Enable = 'off';
                FNum = ScanFileNum(app);
                FNum = FNum + 1;
                NewFileName = sprintf('%02d.mat', FNum);
                FileDuration = TotalTime - LogSt;
                save(fullfile(app.FolderEditField.Value, NewFileName),...
                    'LogRaw','FileDuration','GroundTruth');
%                 app.Label.Text = sprintf('%02d', FNum);
                ScanFileNum(app);
%                 set(handles.CountText,'String', sprintf('%02d', FNum));
                LogRaw = zeros(30,0,3);
                app.DeleteButton.Enable = 'on';
                app.StartButton.Enable = 'on';
            end
        end
%% CSI Ratio Data
        if opts.CSI_Ratio
            DataNext = RawData(:,:,[3 1 2]);
            RatioData = cat(2, RatioData, RawData ./ DataNext);
        end
%% 3.2.1 Multipath Removal
        if opts.MultipathRemoval
            MR_ratio = MRCutTime/0.05335/30;
            % Total bandwith = 20MHz, 64 subcarriers, 30 collected;
            % Assume measured channel = 20MHz * 30 / 64 = 9.375MHz;
            % -9.375~9.375MHz then every interval of time domain is
            % 1/(9.375MHz*2) = 0.05535ms;
            CSI_MR = zeros(size(RawData));
            Start_t = round(size(RawData,1) * MR_ratio);
            for i = 1 : size(RawData,3)
                PDP  = ifft(squeeze(RawData(:,:,i)));
                PDP(Start_t:end,:) = 0;
                CSI_MR(:,:,i) = fft(PDP);
            end
            Mag = cat(2,Mag,abs(CSI_MR));
        else
            Mag = cat(2,Mag,abs(RawData));
        end
        if size(Mag,2) < 200
            continue;
        end
%% Phase calculation
        if opts.PhaseCalculation
            PhaseTmp = zeros(size(RawData));
            for an = 1 : NAntenna
                next = an + 1;
                if next == NAntenna+1
                    next = 1;
                end
                PhaseTmp(:,:,an) = abs(unwrap(angle(...
                    RawData(:,:,next) .* conj(RawData(:,:,an)))));
            end
            Phase1 = cat(2, Phase1, PhaseTmp);
            PhaseTmp = zeros(size(RawData));
            Ph = angle(RawData);
            % Phase calculate - another method for one antenna
            m = cat(2,-28 : 2 : -2, -1, 1, 2 : 2 : 28);
            for an = 1 : NAntenna
                for i = 1 : length(Ph)
                    Tp = zeros(1,NSub);
                    Tp(1) = Ph(1,i,an);
                    flag = 0;
                    for j = 2 : NSub
                        if Ph(j,i,an) - Ph(j-1,i,an) > pi
                            flag = flag + 1;
                        end
                        Tp(j) = Ph(j,i,an) - flag * 2 * pi;
                    end
                    K = (Tp(end) - Tp(1)) / (m(end) - m(1));
                    B = mean(Tp);
                    PhaseTmp(:,i,an) = Tp - K .* m - B; 
                end
            end
            Phase2 = cat(2, Phase2, PhaseTmp);
        end
%% Change the mag size if the total time is lager than 25s
        if TotalTime > MaxPlot + 2 * RefreshInterval
            ChangePlotStartTime = ceil(TotalTime - MaxPlot);
            while ChangePlotStartTime > PlotSt &&...
                    CursorSeries(ChangePlotStartTime) == 1
                ChangePlotStartTime = ChangePlotStartTime - 1;
            end
            ChangeIdx = CursorSeries(ChangePlotStartTime)...
                -CursorSeries(PlotSt);
            Mag(:,1:ChangeIdx,:) = [];
            if opts.PhaseCalculation
                Phase1(:,1:ChangeIdx, :) = [];
                Phase2(:,1:ChangeIdx, :) = [];
            end
            PlotSt = ChangePlotStartTime;
        end
        
        
%% 3.2.2 Hampel Identifier
        if opts.OutlierRemoval
            [~,Idx,~] = hampel...
                (Mag(1,:,HampelAntenna), K0, Nsigma0);
            V = diff([0 Idx 0] == 1);
            St = find(V == 1);
            Ed = find(V == -1);
            TmpIdx = find(Ed - St > 100);
            St = St(TmpIdx);
            Ed = Ed(TmpIdx) - 1;
            OutIdx = [];
            for i = 1 : length(St)
%                 Tmp = Data(sub,St(i):Ed(i),an);
%                 OutlierRM(sub,St(i):Ed(i),an) = Tmp - ...
%                     (mean(Tmp) - Median(St(i):Ed(i)));
                OutIdx = cat(2, OutIdx, St(i):Ed(i));
            end
%             Mag(:,OutIdx,:) = [];
            Outliered = zeros(size(Mag));
            for ham_i = 1:NAntenna
                if ~opts.PCA
                    [Outliered(1,:,ham_i), Hampel_Idx] = hampel(...
                        Mag(1,:,ham_i)',K1,Nsigma1);
                else
                    for ham_j = 1 : NSub
                        Outliered(ham_j,:,ham_i) = hampel(...
                            Mag(ham_j,:,ham_i)',K1,Nsigma1);
                    end
                end
            end
        else
            Outliered = Mag;
        end
        
        Len = size(Mag,2);
        DataDur = TotalTime - PlotSt;
        FS = Len / DataDur;
        if TotalTime < MaxPlot + 2 * RefreshInterval
            inveral = TotalTime /Len;
            T = 0:inveral:(TotalTime - inveral);
        else
            inveral = (TotalTime - PlotSt)/Len;
            T = PlotSt:inveral:(TotalTime - inveral);
        end
%% 3.2.4 PCA
        if opts.PCA % adjust_process >= 3
            PCANeed = zeros(NumPCA,size(Outliered,2),3);
            for an = 1 : NAntenna
                input = squeeze(Outliered(:,:,an));
                input = input';
                coeff = pca(input);
                tmp = input * coeff(:,1:NumPCA);
                PCANeed(:,:,an) = tmp';
            end
        else
            PCANeed = squeeze(Outliered(1:NumPCA,:,:));
%             all(PCA1(:,1)==PCA1(:,3))
        end
%% 3.2.3 Butterworth low-pass filter
        if opts.BandpassFilter
            if abs(FS - 1700) > 500
                FilterFreq = 1700;
            else
                FilterFreq = FS;
            end
             BandPass = fdesign.bandpass('N,F3dB1,F3dB2',...
                FilterOrder,CutoffFreq(1),CutoffFreq(2),FilterFreq);
            % N -- filter order for FIR filters.
            % F3dB1/2 -- cutoff frequency for the point 3 dB point below the 
            % passband value for the first/second cutoff.
            BPFilter = design(BandPass,'butter');
%             Lowpass = fdesign.lowpass('N,F3db',...
%                 FilterOrder, CutoffFreq(2), FS);
%             LPFilter = design(Lowpass,'butter');
%             if TotalTime  < (MaxPlot / 2) ||~exist('LPFilter','var')
%                 LowPass = fdesign.lowpass('N,F3db', FilterOrder,...
%                     CutoffFreq(2), 1700);
%                 LPFilter = design(LowPass,'butter');
%             else
%                 if TotalTime < MaxPlot 
%                     LowPass = fdesign.lowpass('N,F3db',...
%                         FilterOrder, CutoffFreq(2), FS);
%                     LPFilter = design(LowPass,'butter');
%                 end % else use LowPass calculated last time
%             end
%             LPFilter = design(LowPass,'butter');
            Filtered = zeros(NumPCA, Len, NAntenna);
            for i = 1:NAntenna
                if ~opts.PCA
                    Filtered(1,:,i) = filter(BPFilter, PCANeed(1,:,i));
                else
                    for j = 1 : NumPCA
                        Filtered(j,:,i) = filter(BPFilter, PCANeed(j,:,i));
                        Filtered(j,1:ceil(CutTime * FS),i) = ...
                            Filtered(j,ceil(CutTime * FS),i);
                    end
                end
            end
        else
            Filtered = PCANeed;
        end
        
%% 3.2.5 FFT & Action Detection
        if opts.FrequencyDomain
            Fs = FS;
%             T = 1 / Fs;
%             tt = (0:len-1)*T;
            FFTs = fft(csi);
            P2 = abs(FFTs/Len);
            P1 = P2(1:floor(Len/2) + 1);
            P1(2:end-1) = 2 * P1(2:end-1);
            FIdx = Fs * (0:floor(Len/2))/Len;
            WalkFeqIdx = find(FIdx > ActDet(1) & FIdx < ActDet(2));
            FFTEnergy = sum(P1(WalkFeqIdx));
            if TotalTime > ActDet(4) && ~isBlinking &&...
                    FFTEnergy > ActDet(3)
                    %any(P1(WalkFeqIdx)> WalkDet(3))
                isBlinking = true;
                WalkST = TotalTime;
                Slience = false;
            end
            if isBlinking &&...
                    FFTEnergy < ActDet(3)
                    %all(P1(WalkFeqIdx) < WalkDet(3))
                isBlinking = false;
                if BlinkStarted
                    BlinkStarted = false;
                    disp('Walk ended!');
                end
            end
            if ~Slience && isBlinking && ...
                    TotalTime - WalkST > WalkTime(1)
                Slience = true;
                SlienceTime = TotalTime;
                BlinkStarted = true;
%                 disp('Walk started!');
                % Feature Extraction & Classification:
                if HaveClassifier
                    % we do not process the model and recognize people now
%                     [UserD, Confidence] = WiBlink_Detect(Filtered,...
%                         FS, Classifier, ConfidenceThr);
                    if UserD == 0
                        fprintf('Unknown User.\n');
                    else
                        if UserD > length(Users)
                            fprintf('Other user! Ban!!\n');
                        else
                            fprintf('User %s Detected!\n',Users{UserD});
                        end
                    end
                    fprintf('Confidence: %f\n',Confidence);
                else
                    fprintf('.');
                end
            end
            if Slience && isBlinking && ...
                    TotalTime - SlienceTime > WalkTime(1)
                Slience = false;
            end
            
        end
%% Acceleration Estimation
        if opts.AccelerationEstimation
            Data = Filtered;
            SecondDiffs = Data(:,1:DownNum:end,:);
            ThirdDiffs = squeeze(SecondDiffs(1,:,:));
            for an = 1 : NAntenna
                for sub = 1 : NumPCA
                    SecondDiffs(sub, :, an) = SecondDifferentiator...
                        (Data(sub,1:DownNum:end,an),1/FS);
                end
                % twice differentiator
                ThirdDiffs(:,an) = SecondDifferentiator...
                    (SecondDiffs(1,:,an),1/FS);
            end
            for an = 1 : NAntenna
                if opts.PCA
                    [Pks, Idx] = findpeaks(ThirdDiffs(LenSt:LenEd,an), FS,...
                    'MinPeakDistance', MinDist,...
                    'MinPeakHeight', 0.75 * std(Data(LenSt:LenEd,an)));
                else
                    sub = 3;
                    [Pks, Idx] = findpeaks(ThirdDiffs(sub,LenSt:LenEd,an),...
                        FS,'MinPeakDistance', MinDist,...
                        'MinPeakHeight', 0.75 * std(Data(sub,LenSt:LenEd,an)));
                end
            end
        end
%% Mean absolute deviation    
        if opts.MeanAbsoluteDeviation
            Data = squeeze(Filtered(1,:,:));
            MAD = zeros(Len, NAntenna);
            for an = 1 : NAntenna
                for i = 1 + MADwindow/2 : Len - MADwindow/2
                    input = Data(i-MADwindow/2:i+MADwindow/2,an);
                    avg = mean(input);
                    avg_vector = repmat(avg, size(input, 1), size(input, 2));
                    MAD(i,an) = mean(abs(input - avg_vector));
                end
            end
        end
%% Discrete Wavelet Transform
        if opts.DiscreteWaveletTransform
%             DataIn = sum(Filtered(LenSt:LenEd,:), 2);
            Data = squeeze(Filtered(1,LenSt:LenEd,:));
            DWT = zeros(size(Data)); % units: dB
            for an = 1 : NAntenna
                temporary = Data(:,an);
                for k = 1:DWTlevel
                    [Main, temporary] = dwt(temporary, DWTwavelet);
                end
                temporary(1:5) = [];
                temporary(end - 5: end) = [];
                Main(1:5) = [];
                Main(end-5:end) = [];
                PData = 10.*log10(temporary');
%                 PData = abs(PData)./5;
%                 PData = PData - mean(PData);
                DWT(:,an) = PData;
            end
%             STE = energy(PData', 'hamming', 128, 256);
%             STE(length(PData)+1:end) = [];
            DWTinterval = (T(end) - T(1)) / length(temporary);
            DWTidx = T(1) : DWTinterval : T(end) - DWTinterval;
        end
%% 3.3 Draw Figure
        if ~exist('AxesHandle','var')
            AxesHandle = cell(4,1);
        end
        for an = 1 : NAntenna
            ax = subplot(2,2,an, 'Parent', app.Panel);
            AxesHandle{an} = ax;
            plot(ax,T,PCANeed(1,:,an));
            if ~opts.PCA && opts.MultipathRemoval
                plot(ax,T,Mag(1,:,an));
                hold(ax, 'on');
                plot(ax, T,PCANeed(1,:,an));
                plot(ax, T(Hampel_Idx), Mag(1,Hampel_Idx,an), 'sk',...
                    'MarkerSize',0.5);
                hold(ax, 'off');
            end
            if opts.PowerDensitySpectrum
                hold(ax, 'on');
                YLim = get(ax,'ylim');
                [~,w,T,ps] = spectrogram(PCANeed(1,:,an),128,120,128,...
                    FS,'yaxis');
                w = w + YLim(1);
                Drop = find(w > YLim(2));
                w(Drop) = [];
                ps(Drop,:) = [];
                
                h2 = surf(ax,T+DataDur(1),w,10*log10(ps),'edgecolor','none');
                hold(ax, 'off');
                set(ax,'child',[h1 h2]);
                colorbar('East'); % power spectral density
                ax2 = axes('Position',get(ax,'Position'),...
                    'YAxisLocation','right','Color','none');
                YLim2 = [0 YLim(2) - YLim(1)];
                ax2.XAxis.Visible = 'off';
                set(ax2,'xlim', get(ax,'xlim'),'ylim',YLim2);
                ylabel(ax2, 'Frequency/Hz');
                set(gcf,'Child',[ax ax2]);
            end
            if opts.CountAction
                GT = Groundtruth(app,0)+LogSt;
                GT(GT < T(1)) = [];
                axis(ax, 'tight');
                YLim = get(ax, 'ylim');
                hold(ax, 'on');
                for i = 1 : length(GT)
                    plot(ax, [GT(i) GT(i)], YLim, 'k--', 'LineWidth', 2);
                end
                hold(ax, 'off');
            end
            axis(ax, 'tight');
            Xlim = get(ax,'xlim');
            if opts.BandpassFilter
                Xlim(1) = Xlim(1) + (Xlim(2) - Xlim(1)) * 0.1;
                set(ax ,'xlim',Xlim);
            end
%             temp = strfind(file_name,'/');
            temp = DataFileName(1:end);
            temp = strrep(temp,'\','\\');
            temp = strrep(temp,'_','\_');
            if an == NAntenna
                title(ax, ['antenna: ' num2str(an) ' file:' temp]);
            else
                title(ax, ['antenna: ' num2str(an)]);
            end
            xlabel(ax, 'time/s');
            ylabel(ax, 'amplitude/dB');
        end
        if opts.FrequencyDomain
            ax = subplot(2,2,4);
            AxesHandle{4} = ax;
%             Idxs = find(FIdx > WalkDet(1) & FIdx < (WalkDet(2) + 2));
            plot(ax, FIdx(WalkFeqIdx), P1(WalkFeqIdx));
            hold(ax, 'on');
%             line([5 5],[0 FFTEnergy],...
%                 'LineWidth',3);
            hold(ax, 'off');
            set(ax,'YLim',[0 3]);
            xlabel(ax, 'f (Hz)');
            ylabel(ax, '|P1(f)|');
            title(ax, sprintf('Energy: %.3f',FFTEnergy));
        end
        drawnow;
        Duration = toc - TotalTime;
        if (DebugTime)  || (Duration > 1.5 * RefreshInterval)
            display(['duration:' num2str(Duration)]);
        end
        if (Duration <= RefreshInterval)
           pause(RefreshInterval-Duration); 
%            Duration = RefreshInterval;
        end
%         fprintf('start_time: %f\n',plot_start_time);
    end
%     disp(['file was not written for' num2str(refresh_time) 's!...']);
    if opts.LogTotalTime
        temp = strfind(DataFileName,'/');
        log_dir = DataFileName(1:temp(end));
        log_file_name = fullfile(log_dir, 'duration.txt');
        out = fopen(log_file_name,'a');
        fprintf(out,'%s %s\n',log_file_name,num2str(TotalTime));
        fclose(out);
    end
    
%     title(AxesHandle{4},'Ready');
    disp('---------------------------');
    Re = TotalTime;
end
disp('End showing data');
