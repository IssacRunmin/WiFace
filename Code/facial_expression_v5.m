%Facial Expression Data Processing. 
%   Input:
%       Files: is the directory(ies) order you want to process
%       file_choice: if you want to chose to process particullary,select
%       the file, enter the file number. set 0 to process all the files.
%   Output:
%       saved files for each directory
%% Parameters
files = [
 1
%  2 % WiFace_0223_Runmin_Ou_20cm
%  3 % WiFace_0223_Runmin_Ou_30cm
% 15 % WiFace_0319_Runmin_Ou_Office
% 37 % WiFace_0907_Runmin_Ou_90cm
% 38 % WiFace_0907_Runmin_Ou_120cm
% 41 % WiFace_0912_Runmin_Ou_180cm
% 43 % WiFace_0912_Runmin_Ou_60cm
% 45 % WiFace_0914_Runmin_Ou_Cap
% 46 % WiFace_0914_Runmin_Ou_Glasses
% 55 % WiFace_0920_Runmin_Ou_Lap
% 58 % WiFace_0922_Runmin_Ou_Standing
% 80 % WiFace_1107_Runmin_Ou_Back
% 81 % WiFace_1107_Runmin_Ou_Left
];
file_choice = 3;
% Preprocessing
hampel_k = 160;         % Hampel identifier: WindowLen/2
hampel_nsigma = 1.7;    % Hampel identifier: multiplier of sigma
MR_CutTime = 0.5;       % unit: ms (0,1.6005] Multipath Removal Time
cutoff = [0.1 30];      % unit: Hz Bandpass filter cutoff frequency
num_for_PCA = 6;        % reserved number of Principle Component Analysis 
% Segmentation
FitWin = 1024;          % Sliding window size
FitInterval = 32;       % Sliding window step
% threshold_ori = 10;     % Threshold for determining the index of points
length_range = [0.4 3]; % unit: s Discard the waveform that is too long

% information:
description = 'Final-Version';
debug = {               % debug: display file information, draw fig, etc.
    'disp', 1           % display the sampling rate, draw the figure, etc.
    'fig_PCA', 1        % draw the figure with PCA
    'time', 1           % time counting
};         
FeatureNum = 24;
MR_Ratio = MR_CutTime/0.05335/30;  % (500ns, interval, channels)
    % Total bandwith = 20MHz, 64 subcarriers, 30 collected;
    % 802.11n delta_f: Subcarrier freq spacing: 312.5kHz
    % so the measured bandwidth is: 0.3125MHz * 30 = 9.375MHz
    % -9.375~9.375MHz then every interval of time domain is 
    % 1/(9.375MHz*2) = 0.05335ms;
Expressions = {'Happy';
    'Fearful';
    'Surprised';
    'HappilySurprised';
    'AngrilySurprised';
    'FearfullySurprised'};
%% 2. Initialization
file_dir= fullfile('..', 'Data');
dirs = dir(fullfile(file_dir,'WiFace*'));
if length(files) > 1
    file_choice = 0;
end
ExpNum = length(Expressions);
ME_I = 1;
Messages = cell(0,1);
Messages_Idx = cell(0,1);
debug(:,2) = cellfun(@(x) {x==1},debug(:,2));
debug = cell2struct(debug(:,2),debug(:,1),1);
for dir_i = 1 : length(files)
    file_directory = fullfile(file_dir, dirs(files(dir_i)).name);
    file_directory = strcat(file_directory,'/');
    tic
    %% 3. Read file infomation
    % open the file and read the information
    if debug.disp
        disp(['Flie directory: ' file_directory]);
    end
    index = 1;
    file_info = cell(7,30); % maybe no more than 30 files
    file_list = fullfile(file_directory,'file2');
    in = fopen(file_list);
    while ~feof(in)
        tline = fgetl(in);
        content = strsplit(tline, ' ');
        if ~contains(tline,'%') && length(content) > 2
            file_info{1,index} = content{1};% the file name
            file_info{2,index} = str2double(content{2}); % the total time
            file_info{3,index} = str2double(content{3}); % the starting time
            file_info{4,index} = str2double(content{4}); % the sample end time
            if length(content) >= 5
                file_info{5,index} = str2double(content{5}); % unused antenna
                file_info{6,index} = str2double(content{6}); % action count
            else
                file_info{5,index} = str2double(content{5});
            end
            % extract facial expression name / eliminate other signs or numbers
            name_t = file_info{1,index};
            letter_t = isletter(name_t);
            if ~all(letter_t)
                id_t = find(~letter_t);
                name_t = name_t(1:id_t(1)-1);
            end
            idx = strcmpi(name_t, Expressions);
            tmp = find(idx);
            if isempty(tmp)
                warning(['Wrong filename: ' file_info{1,index}]);
                tmp = ExpNum + 1;
            end
            file_info{7,index} = tmp; % file type
            index = index + 1;
        end
    end
    fclose(in);
    file_count = index - 1;
    if file_count < file_choice
        error('no such file choice!');
    end
    % Initialize the directory and cells
    TrainData_All = cell(ExpNum+1,1);
    WaveForms = cell(ExpNum+1, 1);
    endpoint_indices = cell(1, file_count);
    if ~exist(fullfile(file_directory,'csi'), 'dir')
        mkdir(fullfile(file_directory,'csi'));
    end
    if ~exist(fullfile(file_directory,'Figs'), 'dir')
        mkdir(fullfile(file_directory,'Figs'));
    end
    for i = 1 : ExpNum+1
        TrainData_All{i} = zeros(FeatureNum * num_for_PCA,0);
        WaveForms{i} = cell(0,1);
    end
    if (file_choice ~= 0)
        file_info = file_info(:,file_choice);
        file_count = 1;
    else
        file_info = file_info(:,1:file_count);
    end
    clear file_list i in index *_t 
    %% 5. Process the data
    for File = 1:file_count
        try
            file_path = fullfile(file_directory, 'csi', [file_info{1,File}, '.mat']);
            file_information = sprintf('File: %s-%02d-%.3f-%.3f', file_info{1, File},...
                file_info{3, File}, file_info{4, File});
            if debug.disp
                disp(file_information);
            end
            if exist(file_path, 'file')
                load(file_path);
            else
                disp('The first time to read the data...Please wait');
                str = fullfile(file_directory, [file_info{1,File} '.dat']);
                CSI_trace = read_bf_file(str);
                len = size(CSI_trace, 1);
                all_raw_data = zeros(30, len, 3);
                for i = 1:len
                    csi_entry = CSI_trace{i}; %for every packet
                    csi = get_scaled_csi(csi_entry);
                    csi = squeeze(csi(1, :, :)).';
                    all_raw_data(:, i, :) = csi;
                end
                save(file_path,'all_raw_data');
                clear CSI_trace csi_entry csi
            end
            n_antenna = size(all_raw_data,3);
            if debug.time
                fprintf('\tRead data: %.3fs\n',toc);
                tic
            end
            %% Multipath removal
            % Power delay profile
            csi_MR = zeros(size(all_raw_data));
            Start_t = round(size(all_raw_data,1) * MR_Ratio);
            for i = 1 : size(all_raw_data,3)
                PDP = ifft(squeeze(all_raw_data(:,:,i)));
                PDP(Start_t:end, :) = 0;  % make the PDP after Start_t be 0
                csi_MR(:,:,i) = fft(PDP);
            end
            mag = abs(all_raw_data);
            mag_MR = abs(csi_MR);
            if debug.time
                fprintf('\tMultipath removal: %.3fs\n',toc);
                tic
            end
            clear all_raw_data csi_MR Start_t PDP i
            %% Hampel Identifier
            % detect the outliers and replace them with median of neighbor points
            outliered = zeros(size(mag));
            outliered_MR = zeros(size(mag_MR));
            for ham_i = 1:n_antenna
                for ham_j = 1:30
                    outliered(ham_j,:,ham_i) = hampel(mag(ham_j,:,ham_i),...
                        hampel_k,hampel_nsigma);
                    outliered_MR(ham_j,:,ham_i) = hampel(mag_MR(ham_j,:,ham_i),...
                        hampel_k,hampel_nsigma);
                end
            end
            mag = squeeze(mag(1,:,:));
            if debug.time
                fprintf('\tHampel Identifier: %.3fs\n', toc);
                tic
            end
            clear ham_i ham_j mag_MR
            %% Butterworth Filter
            len = size(outliered,2);
            FS = len / (file_info{2,File});
            LowPass = fdesign.lowpass('N,F3db', 10, cutoff(2), FS);
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
            if debug.time
                fprintf('\tButterworth filter: %.3fs\n',toc);
                tic
            end
            if debug.disp
                fprintf('\tSample rate: %.3f\n', sample_rate);
            end
            clear i j FS LowPass LPFilter BandPass BPFilter LowPass2 LPFilter2
            %% PCA
            pca_input = [];
            pca_input_MR = [];% Cutted
            pca1_an = zeros(len, n_antenna);
            pca1_MR_an = zeros(len, n_antenna);
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
            clear pca_input pca_input_MR coeff coeff_MR PCA_input_all
            pca1 = pca1(len_start:len_end);
            outliered = squeeze(outliered(1,:,:));
            outliered_MR = squeeze(outliered_MR(1,:,:));
            filtered = squeeze(filtered(1,:,:));
            filtered_MR = squeeze(filtered_MR(1,:,:));
            if debug.time
                fprintf('\tPCA: %.3fs\n', toc);
                tic
            end
            %% Linear Fitting
            length_range = length_range .* sample_rate;
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
            range = abs(slopes) > 10;
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
                        endpoint_indices{File}(1, j) > length_range(2)/FitInterval
                    endpoint_indices{File}(:,j:j + 1) = [];
                else
                    j = j + 2;
                end
            end
            if length(endpoint_indices{File}) > 3
                temp = find(endpoint_indices{File}(1,2:2:end) - endpoint_indices{File}...
                    (1,1:2:end) < length_range(1)/FitInterval);
                endpoint_indices{File}(:,[temp * 2 - 1 , temp * 2]) = [];
            end
            if ~isempty(endpoint_indices{File})
                endpoint_indices{File}(1,:) = Idx(endpoint_indices{File}(1,:));
            end

            if debug.disp
                fprintf('\tFeature number: %d\n', size(endpoint_indices{File},2) / 2);
            end
            if debug.time
                fprintf('\tLinear fitting: %.3fs\n',toc);
                tic
            end

            clear Index j p range s I temp
            %% Feature Extraction
            ActNum = floor(size(endpoint_indices{File}, 2) / 2);
            if ActNum == 0 
                warning('No Segmentation');
            else
                Ends = reshape(endpoint_indices{File}(1,:), 2, ActNum);
                Ends(1,:) = floor(Ends(1,:) - 0.4 * sample_rate);
                Ends(2,:) = ceil(Ends(2,:) + 0.4 * sample_rate);
                if Ends(1,1) < 1
                    Ends(1,1) = 1;
                end
                if Ends(2,end) > length(PCA_need)
                    Ends(2,end) = length(PCA_need);
                end
                % Feature Exaction
                TrainData = zeros(FeatureNum * num_for_PCA, ActNum);
                for PCA_i = 1 : num_for_PCA
                    Feature_t = NaN(FeatureNum, ActNum);
                    for j = 1 : ActNum
                        Wave = PCA_need(Ends(1,j):Ends(2,j), PCA_i);
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
                        Feature_t(1,j) = (Wave(1) + Wave(end)) / 2; % avg of start and end
                        Feature_t(2,j) = time(Ends(2,j))-time(Ends(1,j)); % time
                        Feature_t(3,j) = yf(1); % power of time domain
                        Feature_t(4,j) = 0; % first peak of feq
                        Feature_t(5,j) = 0; % second peak of feq
                        if length(pksf) > 1
                            Feature_t(4,j) = pksf(1);
                            Feature_t(5,j) = pksf(2);
                        else
                            if length(pksf) == 1
                                Feature_t(4,j) = pksf(1);
                            end
                        end
                        Feature_t(6,j) = mean(Wave); % avg of all points;
                        Feature_t(7,j) = max(Wave); % max
                        Feature_t(8,j) = min(Wave); % min
                        Feature_t(9,j) = prctile(Wave,90); % max90th
                        Feature_t(10,j) = prctile(Wave,10); % min10th
                        Feature_t(11,j) = Feature_t(7,j) - Feature_t(8,j); % Range
                        Feature_t(12,j) = std(Wave); % standard deviation
                        Feature_t(13,j) = skewness(Wave); % skewness
                        Feature_t(14,j) = kurtosis(Wave); % kurtosis
            %             TrainData_t(15,j) = getEntropy(Wave);
                        Feature_t(15,j) = 0;
                        Feature_t(16,j) = Feature_t(12,j) / Feature_t(6,j) * 1000;
                        % CV: std/mean * 1000
                        Feature_t(17,j) = median(Wave);% Median
                        Feature_t(18,j) = prctile(Wave,75);% Q3
                        Feature_t(19,j) = prctile(Wave,25);% Q1
                        Feature_t(20,j) = Feature_t(18,j)-Feature_t(19,j); % Q3-Q1
                        Feature_t(21,j) = [0.25 0.5 0.25] * ...
                            prctile(Wave,[25 50 75])'; % SM(...Mean)
                        Feature_t(22,j) = mean(Wave .^ 2); % E(X^2)
                        Feature_t(23,j) = CovMatrix(1,2); % Cov(act, time)
                        Feature_t(24,j) = CorrCoefMatrix(1,2); 
                    end
                    TrainData((PCA_i-1)*FeatureNum+1:PCA_i*FeatureNum,:)=Feature_t;
                end
                if ~isempty(TrainData)
                    TrainData(:,1) = [];
                end
                TrainData_All{file_info{7,File}} = cat(2,...
                    TrainData_All{file_info{7,File}}, TrainData);

                clear ActNum TrainData_t tmp TrainData PCA_i Ends j Wave N n t yf 
                clear magn f pksf Y_EY CovMatrix CorrCoefMatrix 
            end
            if debug.time
                fprintf('\tFeature extraction: %.3fs\n',toc);
                tic
            end
            draw_figure_v5;
        catch ME
            Messages{ME_I,1} = getReport(ME);
            Messages_Idx{ME_I,1} = fullfile(file_directory, file_information);
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
    FeatureData = TrainData_All;
    clear tmp all_raw_data csiall File
    clear FS interval len
    clear time TrainData_All mag_MR outliered_MR
    clear filtered_MR u_antenna n_antenna time_series
    SaveTime = datestr(datetime(),'mmddHHMMSS');
    if file_choice == 0
        save([file_directory '/TrainData_' description '_' SaveTime '.mat'],...
            'FeatureData');
    end
end