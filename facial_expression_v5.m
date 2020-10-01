%   input:
%       file_directory is the data (set) you want to process
%       file_choice: if you want to chose to process particullary,select
%       the file,enter the file number. No this parmeter to process all.
%   output:
%       saved files for each directory
%% 1. Initialize variables
save_endpoint = true;
endpoint_loaded = false;
Discription = 'BandPass-AllCutoff';
Files = [
 2 % WiFace_0223_Runmin_Ou_20cm
 3 % WiFace_0223_Runmin_Ou_30cm
15 % WiFace_0319_Runmin_Ou_Office
37 % WiFace_0907_Runmin_Ou_90cm
38 % WiFace_0907_Runmin_Ou_120cm
41 % WiFace_0912_Runmin_Ou_180cm
42 % WiFace_0912_Runmin_Ou_240cm
43 % WiFace_0912_Runmin_Ou_60cm
45 % WiFace_0914_Runmin_Ou_Cap
46 % WiFace_0914_Runmin_Ou_Glasses
55 % WiFace_0920_Runmin_Ou_Lap
58 % WiFace_0922_Runmin_Ou_Standing
80 % WiFace_1107_Runmin_Ou_Back
81 % WiFace_1107_Runmin_Ou_Left
];
file_choice = 1; % 8 
if length(Files) > 1
    file_choice = 0;
end
if file_choice ~= 0
    endpoint_loaded = false;
    save_endpoint = false;
end

% debug
debug = 1; % display the sampling rate, draw the figure, ect.
debug_PCA = 1; % draw the figure with PCA
debug_time = 1; % time counting
hampel_k = 160;%200;
hampel_nsigma = 1.7;
MR_CutTime = 0.5;       % Multipath Removal Time
% 0 < MR_CutTime < 1.6005
cutoff = [0.1 30];
cutoff_all = 10:10:80;
cutoff_low = 20;
max_samples = 80;
% linear fitting
FitWin = 1024; % window size
FitInterval = 32;
threshold_ori = 10; % Threshold for determining the index of points
biggest_length = 8192; % 4096; % discard the waveform that is too long
smallest_length = 0.4; % s
extend = [0.4 0.4]; % s, segmentation point before & after expression
% PCA
    num_for_PCA = 20;
% DWT
DWT_level = 0;
DWT_wavelet = 'db1';

% file information:
% NFold = 10;
FilterLen = length(cutoff_all);
max_file = 30;
ExpNum = 6;
FeatureNum = 24;
num_subcarrier = 30;
MR_Ratio = MR_CutTime/0.05335/30;  % (500ns, interval, channels)
    % Total bandwith = 20MHz, 64 subcarriers, 30 collected;
    % 802.11n delta_f: Subcarrier freq spacing: 312.5kHz
    % so the measured bandwidth is: 0.3125MHz * 30 = 9.375MHz
    % -9.375~9.375MHz then every interval of time domain is 
    % 1/(9.375MHz*2) = 0.05335ms;
Expressions = {'happy';
    'fearful';
    'surprised';
    'happilysurprised';
    'angrilysurprised';
    'fearfullysurprised'};
SegFile = 'Segmentation_v3.mat';
ME_I = 1;
Messages = cell(0,1);
Messages_Idx = cell(0,1);
%% 2. Initialization
FileDir= 'WiFace_Data';
Dirs = dir(fullfile(FileDir,'WiFace*'));
for Dir_i = 1 : length(Files)


file_directory = fullfile(FileDir, Dirs(Files(Dir_i)).name);
file_directory = strcat(file_directory,'/');
tic
%% 3. Read file infomation
% open the file and read the information
if debug == 1
    disp(['Flie directory: ' file_directory]);
end
index = 1;
file_info = cell(7,max_file);
file_list = strcat(file_directory,'file2');
in = fopen(file_list);
while ~feof(in)
    tline = fgetl(in);
    if ~contains(tline,'%') && length(strfind(tline,' ')) > 2
        interval = strfind(tline,' ');
        file_info{1,index} = tline(1:interval(1)-1);% the file name
        file_info{2,index} = str2double(tline(interval(1):interval(2)));% the total time
        file_info{3,index} = str2double(tline(interval(2):interval(3)));% the sample starting time
        file_info{4,index} = str2double(tline(interval(3):interval(4)));% the sample end time
        if length(interval) >= 5
            file_info{5,index} = str2double(tline(interval(4):interval(5)));% unused antenna
            file_info{6,index} = str2double(tline(interval(5):end));% action times
        else
            file_info{5,index} = str2double(tline(interval(4):end));% unused antenna
        end
        name_t = file_info{1,index};
        letter_t = isletter(name_t);
        if ~all(letter_t)
            id_t = find(~letter_t);
            name_t = name_t(1:id_t(1)-1);
        end
        name_t = lower(name_t);
        switch name_t
            case Expressions{1}
                tmp = 1;
            case Expressions{2}
                tmp = 2;
            case Expressions{3}
                tmp = 3;
            case Expressions{4}
                tmp = 4;
            case Expressions{5}
                tmp = 5;
            case Expressions{6}
                tmp = 6;
            otherwise
                warning(['Wrong filename: ' file_info{1,index}]);
                tmp = 7;
        end
        file_info{7,index} = tmp; % file type
%         file_info{5,index} = str2double(tline(interval(4):length(tline)));% the antenna unused
        index = index + 1;
    end
    if (index == max_file + 1)
        break;
    end
end
fclose(in);
file_count = index - 1;

%Process file information
% if I just choose one file,just process on file and do not process KNN
if file_count < file_choice
    error('no such file choice!');
end
TrainData_All = cell(ExpNum+1,1);
TrainData_Cutoff = cell(FilterLen,ExpNum+1);
WaveForms = cell(ExpNum+1, 1);
endpoint_indices = cell(1, file_count);
if ~exist([file_directory,'csi'], 'dir')
    mkdir([file_directory,'csi']);
end
if ~exist([file_directory,'Figs'], 'dir')
    mkdir([file_directory,'Figs']);
end
for f = 1 : length(cutoff_all)
    for i = 1 : ExpNum+1
        TrainData_All{i} = zeros(FeatureNum * num_for_PCA,0);
        WaveForms{i} = cell(0,1);
        TrainData_Cutoff{f,i} = zeros(FeatureNum * num_for_PCA,0);
    end
%     TrainData_Cutoff{f} = TrainData_All;
end
if (file_choice ~= 0)
    file_info = file_info(:,file_choice);
    file_count = 1;
else
    file_info = file_info(:,1:file_count);
end

% if debug_time == 1
%     disp(['    File information:' num2str(toc)]);
%     tic
% end

clear file_list i in index *_t 
%% 5. Process the data
for File = 1:file_count
try
    str = strcat([file_directory,'csi/'],[ file_info{1,File},'.mat']);
    SegFileName = [file_info{1, File} '-' num2str(file_info{3, File}) ...
        '-' num2str(file_info{4, File})];
%     SegFileName = file_info{1, File};
    if debug == 1
        disp(['File: ' SegFileName]);
    end
    try
%         error(' ');
        load(str);
        
    catch
        disp('The first time to read the data...Please wait');
        str = [file_directory  file_info{1,File} '.dat'];
        CSI_trace = read_bf_file(str);
        len = size(CSI_trace, 1);
        all_raw_data = zeros(num_subcarrier, len, 3);
        for i = 1:len
            csi_entry = CSI_trace{i}; %for every packet
            csi = get_scaled_csi(csi_entry);
            csi = squeeze(csi(1, :, :)).';
            all_raw_data(:, i, :) = csi;
        end
        str = strcat([file_directory,'csi/'],[ file_info{1,File},'.mat']);
        save(str,'all_raw_data');
        clear CSI_trace csi_entry csi
    end
    n_antenna = size(all_raw_data,3);
    if debug_time == 1
        disp(['    Read Data:' num2str(toc)]);
        tic
    end
    Parameters = cat(2, hampel_k, hampel_nsigma, cutoff_low, ...
    FitWin,u_antenna);
    %% 5.1 Multipath removal
    % Power delay profile
    ProcessData = cell(length(cutoff_all),1);
    csi_MR = zeros(size(all_raw_data));
    Start_t = round(size(all_raw_data,1) * MR_Ratio);
    for i = 1 : size(all_raw_data,3)
        PDP = ifft(squeeze(all_raw_data(:,:,i)));
        PDP(Start_t:end, :) = 0;
        csi_MR(:,:,i) = fft(PDP);
    end
    mag = abs(all_raw_data);
    mag_MR = abs(csi_MR);
    if debug_time == 1
        disp(['    Multipath Removal:' num2str(toc)]);
        tic
    end
    clear all_raw_data csi_MR Start_t PDP i
    %% 5.2 Hampel Identifier
    outliered = zeros(size(mag));
    outliered_MR = zeros(size(mag_MR));
    for ham_i = 1:n_antenna
        for ham_j = 1:num_subcarrier
            outliered(ham_j,:,ham_i) = hampel(mag(ham_j,:,ham_i),...
                hampel_k,hampel_nsigma);
            outliered_MR(ham_j,:,ham_i) = hampel(mag_MR(ham_j,:,ham_i),...
                hampel_k,hampel_nsigma);
        end
    end
    mag = squeeze(mag(1,:,:));
    mag_MR = squeeze(mag_MR(1,:,:));
    if debug_time == 1
        disp(['    Hampel Identifier:' num2str(toc)]);
        tic
    end
    clear ham_i ham_j
    %% 5.3 Butterworth Filter-- bandpass filter
    len = size(outliered,2);
    FS = len / (file_info{2,File});
    LowPass = fdesign.lowpass('N,F3db', 10, cutoff_low, FS);
    LPFilter = design(LowPass, 'butter');
    BandPass = fdesign.bandpass('N,F3dB1,F3dB2',10,cutoff(1),cutoff(2),FS);
    BPFilter2 = design(BandPass,'butter');
        % N -- filter order for FIR filters.
        % F3dB1 -- cutoff frequency for the point 3 dB
        % point below the passband value for the first cutoff.
        % F3dB2 -- cutoff frequency for the point 3 dB 
        % point below the passband value for the second cutoff.
    filtered = zeros(size(outliered));
    filtered_MR = zeros(size(outliered_MR));
    for i = 1:size(outliered, 3)
        for j = 1:size(outliered, 1)
            filtered(j,:,i) = filter(LPFilter, outliered(j,:,i));
            filtered_MR(j,:,i) = filter(BPFilter2, outliered_MR(j,:,i));
        end
    end
    
    sample_rate = len / file_info{2,File};
    plot_interval = 1 / sample_rate;
    len_t = length(filtered);
    FS_t = len_t / (file_info{4,File} - file_info{3,File});
    Interval = 1 / FS_t;
    time = file_info{3,File}:Interval:(file_info{4,File} - Interval);
    len_start = floor(len * file_info{3,File} / file_info{2,File}) + 1;
    len_end = min(floor(len * file_info{4,File} / file_info{2,File}), len);
    if debug_time == 1
        disp(['    Butterworth filter:' num2str(toc)]);
        tic
    end
    if debug == 1
        disp(['    Sample rate: ' num2str(sample_rate)]);
    end
    clear i j FS LowPass LPFilter BandPass BPFilter LowPass2 LPFilter2
    %% 5.4 PCA
    pca_input = [];
    pca_input_MR = [];% Cutted
    pca1_an = zeros(len, n_antenna);
    pca1_MR_an = zeros(len, n_antenna);
    PCA_input_all = cell(FilterLen,1);
    for j = 1:n_antenna
        pca_input = cat(1, pca_input, squeeze(filtered(:,:,j)));
        pca_input_MR = cat(1, pca_input_MR,...
            squeeze(filtered_MR(:,len_start:len_end,j)));
        tmp = squeeze(filtered(:,:,j))';
        coeff = pca(tmp);
        pca1_an(:,j) = tmp * coeff(:,1);
        tmp = squeeze(filtered_MR(:,:,j))';
        coeff = pca(tmp);
        pca1_MR_an(:,j) = tmp * coeff(:,1);
    end
    pca_input = pca_input.';
    pca_input_MR = pca_input_MR.';
    coeff = pca(pca_input);
    pca1 = pca_input * coeff(:,1);
    coeff_MR = pca(pca_input_MR);
    PCA_need = pca_input_MR * coeff_MR(:,1:num_for_PCA);
    
    for f = 1 : FilterLen
        BandPass = fdesign.bandpass('N,F3dB1,F3dB2',10,...
                    cutoff(1),cutoff_all(f),sample_rate);
        BPFilter2 = design(BandPass,'butter');
        filtered_tmp = zeros(size(filtered));
        for i = 1:size(filtered, 3)
            for j = 1:size(filtered, 1)
                filtered_tmp(j,:,i) = filter(BPFilter2, outliered_MR(j,:,i));
            end
        end
        pca_input_cutoff = [];
        for j = 1:n_antenna
            pca_input_cutoff = cat(1, pca_input_cutoff,...
                squeeze(filtered_tmp(:,len_start:len_end,j)));
        end
        pca_input_cutoff = pca_input_cutoff';
        coeff = pca(pca_input_cutoff);
        ProcessData{f} = pca_input_cutoff * coeff(:,1:num_for_PCA);
    end
    clear pca_input pca_input_MR coeff coeff_MR PCA_input_all
    pca1 = pca1(len_start:len_end);
    outliered = squeeze(outliered(1,:,:));
    outliered_MR = squeeze(outliered_MR(1,:,:));
    filtered = squeeze(filtered(1,:,:));
    filtered_MR = squeeze(filtered_MR(1,:,:));
    if debug_time == 1
        disp(['    PCA:' num2str(toc)]);
        tic
    end
%     coreef = corrcoef(pca_input_MR);
%     s = spectrogram(pca1);
    %% 5.5 Linear Fitting
    try 
        if ~endpoint_loaded
            error('linear fitting');
        end
        load([file_directory SegFile]);
    catch
        Idx = FitWin / 2:FitInterval:length(pca1) - FitWin / 2;
        slopes = zeros(1,length(Idx));
        temp = 10000 * zscore(pca1);
        if iscolumn(temp)
            temp = temp.';
        end
        for j = 1 : length(Idx)
            Idx_tt = Idx(j) - FitWin / 2 + 1:Idx(j) + FitWin / 2;
            p = polyfit(Idx_tt, temp(Idx_tt),1);
            slopes(j) = p(1);    
        end
        endpoint_indices{File} = [];
        range = abs(slopes) > threshold;
        range(1) = 0;
        range(end) = 0;
        for j = 2:length(range)
            if range(j) == 1 && range(j - 1) == 0
                s = j;
            end
            if range(j) == 0 && range(j - 1) == 1
                t = j - 1;
                [~,I] = max(abs(slopes(s:t)));
                I = I + s - 1;
                if slopes(I) > 0
                    endpoint_indices{File} = cat(2, endpoint_indices{File}, [I;1]);
                else
                    endpoint_indices{File} = cat(2, endpoint_indices{File}, [I;-1]);
                end
            end
        end
        j = 2;
        while j <= size(endpoint_indices{File}, 2)
            if endpoint_indices{File}(2, j) == endpoint_indices{File}(2, j - 1)
                if endpoint_indices{File}(1, j) - ...
                        endpoint_indices{File}(1, j - 1) < 1800/FitInterval
                    endpoint_indices{File}(:,j - 1) = [];
                else
                    endpoint_indices{File}(:,j) = [];
                end
            else
                j = j + 1;
            end
        end
        if rem(length(endpoint_indices{File}), 2) == 1
            temp = length(endpoint_indices{File});
            endpoint_indices{File}(:, temp) = [];
        end
        j = 1;
        while j < size(endpoint_indices{File}, 2)
            if endpoint_indices{File}(1, j + 1) - ...
                    endpoint_indices{File}(1, j) > biggest_length/FitInterval
                endpoint_indices{File}(:,j:j + 1) = [];
            else
                j = j + 2;
            end
        end
        if length(endpoint_indices{File}) > 3
            temp = find(endpoint_indices{File}(1,2:2:end) - endpoint_indices{File}...
                (1,1:2:end) < smallest_length * sample_rate/FitInterval);
            endpoint_indices{File}(:,[temp * 2 - 1 , temp * 2]) = [];
        end
        if size(endpoint_indices{File},2) > max_samples * 2 && file_choice == 0
            endpoint_indices{File}(:,max_samples * 2 +...
                1:size(endpoint_indices{File},2)) = [];
        end
        Tmp_RecordPeakDur;
        if ~isempty(endpoint_indices{File})
            endpoint_indices{File}(1,:) = Idx(endpoint_indices{File}(1,:));
        end
        if file_choice == 0 && save_endpoint
            save([file_directory SegFile],'endpoint_indices');
        end
    end
    if debug == 1
        disp(['    Feature number: ' num2str(size(endpoint_indices{File},2) / 2)]);
    end
    if debug_time == 1
        disp(['    Linear fitting:' num2str(toc)]);
        tic
    end
    
    clear Index j p range s I temp
    %% 5.5 Feature Extraction
    
    ActNum = floor(size(endpoint_indices{File}, 2) / 2);
    if ActNum == 0 
        warning('No Segmentation');
    else
        Ends = reshape(endpoint_indices{File}(1,:), 2, ActNum);
        Ends(1,:) = floor(Ends(1,:) - extend(1) * sample_rate);
        Ends(2,:) = ceil(Ends(2,:) + extend(2) * sample_rate);
        if Ends(1,1) < 1
            Ends(1,1) = 1;
        end
        if Ends(2,end) > length(PCA_need)
            Ends(2,end) = length(PCA_need);
        end
        % DWT 
        WaveForm_t = cell(ActNum, 1);
        for j = 1 : ActNum
            tmp = NaN(num_for_PCA, ceil(biggest_length / (2^DWT_level)));
            for PCA_i = 1 : num_for_PCA
                PCA_t = PCA_need(Ends(1,j):Ends(2,j), PCA_i);
                if isrow(PCA_t)
                    PCA_t = PCA_t';
                end
                for k = 1 : DWT_level
                    [PCA_t, ~] = dwt(PCA_t, DWT_wavelet, 'mode', 'per');
                end
                tmp(PCA_i,1:length(PCA_t)) = PCA_t;
            end
            tmp(:, end-sum(isnan(tmp(1,:)))+1 : end) = [];
            WaveForm_t{j} = tmp;
        end
        Exp_i = file_info{7,File};
        WaveForms{Exp_i} = cat(1, WaveForms{Exp_i}, WaveForm_t);
        % Feature Exaction
        for Fi = 1 : FilterLen
            TrainData = NaN(FeatureNum * num_for_PCA, ActNum);
            for PCA_i = 1 : num_for_PCA
                TrainData_t = NaN(FeatureNum, ActNum);

                for j = 1 : ActNum
                    Wave = ProcessData{Fi}(Ends(1,j):Ends(2,j), PCA_i);
                    % FFT
                    N = length(Wave);
                    n = 0 : N-1;
                    t = n / sample_rate;
                    yf = fft(Wave,N);
                    magn = abs(yf(1:ceil(length(yf)/2)));
                    f = n .* sample_rate ./ N;
                    [pksf,~] = findpeaks(magn, 'MINPEAKHEIGHT', 0);
                    Y_EY = time(Ends(1,j):Ends(2,j)) - mean(time(Ends(1,j):Ends(2,j)));
                    CovMatrix = cov(Wave,Y_EY);
                    CorrCoefMatrix = corrcoef(Wave,Y_EY);
                    TrainData_t(1,j) = (Wave(1) + Wave(end)) / 2; % avg of start and end
                    TrainData_t(2,j) = time(Ends(2,j))-time(Ends(1,j)); % time
                    TrainData_t(3,j) = yf(1); % power of time domain
                    TrainData_t(4,j) = 0; % first peak of feq
                    TrainData_t(5,j) = 0; % second peak of feq
                    if length(pksf) > 1
                        TrainData_t(4,j) = pksf(1);
                        TrainData_t(5,j) = pksf(2);
                    else
                        if length(pksf) == 1
                            TrainData_t(4,j) = pksf(1);
                        end
                    end
                    TrainData_t(6,j) = mean(Wave); % avg of all points;
                    TrainData_t(7,j) = max(Wave); % max
                    TrainData_t(8,j) = min(Wave); % min
                    TrainData_t(9,j) = prctile(Wave,90); % max90th
                    TrainData_t(10,j) = prctile(Wave,10); % min10th
                    TrainData_t(11,j) = TrainData_t(7,j) - TrainData_t(8,j); % Range
                    TrainData_t(12,j) = std(Wave); % standard deviation
                    TrainData_t(13,j) = skewness(Wave); % skewness
                    TrainData_t(14,j) = kurtosis(Wave); % kurtosis
        %             TrainData_t(15,j) = getEntropy(Wave);
                    TrainData_t(15,j) = 0;
                    TrainData_t(16,j) = TrainData_t(12,j) / TrainData_t(6,j) * 1000;
                    % CV: std/mean * 1000
                    TrainData_t(17,j) = median(Wave);% Median
                    TrainData_t(18,j) = prctile(Wave,75);% Q3
                    TrainData_t(19,j) = prctile(Wave,25);% Q1
                    TrainData_t(20,j) = TrainData_t(18,j)-TrainData_t(19,j); % Q3-Q1
                    TrainData_t(21,j) = [0.25 0.5 0.25] * ...
                        prctile(Wave,[25 50 75])'; % SM(...Mean)
                    TrainData_t(22,j) = mean(Wave .^ 2); % E(X^2)
                    TrainData_t(23,j) = CovMatrix(1,2); % Cov(act, time)
                    TrainData_t(24,j) = CorrCoefMatrix(1,2); 
                end
                TrainData((PCA_i-1) * FeatureNum + 1 : PCA_i * FeatureNum,:) = TrainData_t;
            end
            if ~isempty(TrainData)
                TrainData(:,1) = [];
            end
            tmp = TrainData_Cutoff{Fi, file_info{7,File}};
        %     tmp = cat(2,tmp,zscore(TrainData));
            tmp = cat(2,tmp,TrainData);
            TrainData_Cutoff{Fi, file_info{7,File}} = tmp;
        end
    end
    clear ActNum TrainData_t tmp TrainData PCA_i Ends j Wave N n t yf 
    clear magn f pksf Y_EY CovMatrix CorrCoefMatrix 
catch ME
    Messages{ME_I,1} = getReport(ME);
    Messages_Idx{ME_I,1} = fullfile(file_directory, SegFileName);
    if file_choice ~= 0
        rethrow(ME);
    end
    fprintf('\n\n\n___________ERROR_START__________\n');
    disp(Messages_Idx{ME_I});
    disp(Messages{ME_I});
    fprintf('___________ERROR_END____________\n\n\n');
    ME_I = ME_I + 1;
end
end
TrainData = TrainData_All;
clear tmp all_raw_data csiall File
clear FS interval len
clear time TrainData_All mag_MR outliered_MR
clear filtered_MR u_antenna n_antenna time_series
SaveTime = datestr(datetime(),'mmddHHMMSS');
save([file_directory 'TrainData_' Discription '_' SaveTime '.mat'],'TrainData_Cutoff',...
    'Discription','WaveForms');
end