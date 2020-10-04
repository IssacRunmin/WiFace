% Options
BalanceDataset = false;
NFold_All = 10;
Once = false;
RandomDataPrem = false;
MapMinMax = false;

PCAs = [1 2 3 4 5 8 11 14 17 20];
num_for_PCA = 3;
cutoff_all = 10:10:80;
FeatureNum = 24;
Files = [
% 2 % WiFace_0223_Runmin_Ou_20cm
%  3 % WiFace_0223_Runmin_Ou_30cm
% 15 % WiFace_0319_Runmin_Ou_Office
37 % WiFace_0907_Runmin_Ou_90cm
% 38 % WiFace_0907_Runmin_Ou_120cm
% 41 % WiFace_0912_Runmin_Ou_180cm
% 42 % WiFace_0912_Runmin_Ou_240cm
% 43 % WiFace_0912_Runmin_Ou_60cm
% 45 % WiFace_0914_Runmin_Ou_Cap
% 46 % WiFace_0914_Runmin_Ou_Glasses
% 55 % WiFace_0920_Runmin_Ou_Lap
% 58 % WiFace_0922_Runmin_Ou_Standing
% 80 % WiFace_1107_Runmin_Ou_Back
% 81 % WiFace_1107_Runmin_Ou_Left
];
Users = [4 6 8 9 10 12 13 14 16 26 27 37 51 57 65 70 71 73 78];
ULen = length(Users);
PLen = length(PCAs);
TrainTime = zeros(ULen, NFold_All, PLen);
TrainCPU = zeros(ULen, NFold_All, PLen);
TestTime = zeros(ULen, NFold_All, PLen);
TestCPU = zeros(ULen, NFold_All, PLen);
FileDir = 'WiFace_Data';
ME_I = 1;
Messages = cell(1);
Messages_Idx = cell(1);
Expressions = {'happy';
    'fearful';
    'surprised';
    'happilysurprised';
    'angrilysurprised';
    'fearfullysurprised'};
ExpNum = length(Expressions);
Dirs = dir(fullfile(FileDir,'WiFace*'));

for Dir_i = 1 : length(Files)
    try
    Dir = dir(fullfile(FileDir, Dirs(Files(Dir_i)).name,...
        'TrainData_BandPass-AllCutoff*.mat'));
    load(fullfile(Dir(end).folder,Dir(end).name), 'TrainData_Cutoff');
    TrainData_Cutoff(:,ExpNum+1) = [];
    if MapMinMax
        for i = 1 : size(TrainData_Cutoff,1)
            for j = 1 : size(TrainData_Cutoff,2)
                TrainData_Cutoff{i,j} = mapminmax(TrainData_Cutoff{i,j});
                % For each row(all expressions), map to -1 ~ 1
            end
        end
    end
%     for i = 1 : size(TrainData_Cutoff,1)
%         for j = 1 : size(TrainData_Cutoff,2)
%             TrainData_Cutoff{i,j}(num_for_PCA * FeatureNum + 1 : end,:) = [];
%             % For each row(all expressions), map to -1 ~ 1
%         end
%     end
    TrainData = TrainData_Cutoff(3,:)';
    ActNums = zeros(length(TrainData),1);
    NFold = NFold_All;
    for i = 1 : length(TrainData)
        ActNums(i) = size(TrainData{i},2);
    end
    if length(ActNums) > ExpNum
        ActNums(ExpNum + 1 : end) = [];
        TrainData(ExpNum + 1 : end) = [];
    end
    if BalanceDataset
        ActNums(ActNums == 0) = [];
        Min = min(ActNums);
        for i = 1 : length(TrainData)
            if ~isempty(TrainData{i}) && length(TrainData{i}) > Min
                TrainData{i}(:,Min+1:end) = [];
            end
        end
    else
        for t = 1 : ExpNum
            if size(TrainData{t},2) > 50
                TrainData{t}(:,51:end) = [];
            end
        end
    end
    
    Num_Sample_group = zeros(NFold_All, ExpNum);
    TrainDataLabel = cell(ExpNum,1);
    for t = 1 : ExpNum
        TrainDataLabel{t} = ones(1,size(TrainData{t},2));
        yushu = rem(length(TrainDataLabel{t}), NFold);
        shang = (length(TrainDataLabel{t}) - yushu) / NFold;
        Num_Sample_group(:, t) = ones(NFold, 1) * shang;
        Num_Sample_group(1:yushu, t) = Num_Sample_group(1:yushu, t) + 1;
        if RandomDataPrem
            perm = randperm(size(TrainData{t},2));
            TrainData{t} = TrainData{t}(:,perm);
        end
        
    end
    TrainData_All_PCA = TrainData;
    % Ten fold cross validation
    if Once
        NFold = 1;
    end
    All_Accuracy = cell(length(num_for_PCA),1);
    All_Confusion = cell(length(num_for_PCA),1);
    All_Template = cell(length(num_for_PCA),1);
    SaveTime = datestr(datetime(),'mmddHHMMSS');
    for PCA_i = 1 : length(num_for_PCA)
        TrainData_All = TrainData_All_PCA;
        for i = 1 : length(TrainData_All)
            TrainData_All{i}(num_for_PCA(PCA_i) * FeatureNum + 1 : end,:) = [];
        end
        RoundAccuracy = zeros(1,NFold);
        confusion_matrices = zeros(ExpNum, ExpNum, NFold);
        start_indices = ones(1, ExpNum);
        PredictLabelTotal = cell(ExpNum, 1);
        for Round = 1 : NFold
            TrainData = TrainData_All;
            TrainLabel = TrainDataLabel;
            TestData = cell(ExpNum,1);
            % Prepare the train and test sample:
            x_test = [];
            y_test = [];
            x_train = [];
            y_train = [];
            for EventCount = 1 : ExpNum
                samples = TrainData{EventCount}';
                s = start_indices(EventCount);
                t = s + Num_Sample_group(Round, EventCount) - 1;
                x_test = [x_test;samples(s:t,:)];
                y_test = [y_test;EventCount * ones(t - s + 1,1)];
                samples(s:t,:) = [];
                x_train= [x_train;samples];
                y_train = [y_train;EventCount * ones(size(samples, 1),1)];
                start_indices(EventCount) = t + 1;
            end
            CPUSt = cputime();
            tic
            if Round == 1
                ClassifierTemplate = fitcensemble(x_train,y_train,'OptimizeHyperparameters','auto',...
                    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
                    'expected-improvement-plus','MaxObjectiveEvaluations',20,...
                    'Kfold',5,'Verbose',1,'ShowPlots',false, 'UseParallel', true));
                Learner = ClassifierTemplate.ModelParameters.LearnerTemplates{1};
            end
            
            if strcmp(ClassifierTemplate.Method, 'Bag')
                Classifier = fitcensemble(x_train, y_train, 'Method',...
                    ClassifierTemplate.Method, 'Learners',Learner,...
                    'NumLearningCycles',ClassifierTemplate.NumTrained);
            else
                Classifier = fitcensemble(x_train, y_train, 'Method',...
                    ClassifierTemplate.Method, 'Learners',Learner,...
                    'NumLearningCycles',ClassifierTemplate.NumTrained,...
                    'LearnRate',ClassifierTemplate.ModelParameters.LearnRate);
            end
            TrainTime(Dir_i, Round, PCA_i) = toc;
            TrainCPU(Dir_i, Round, PCA_i) = cputime() - CPUSt;
        %     kflc = kfoldLoss(Classifier,'Mode','cumulative');
        %     figure;
        %     plot(kflc);
        %     ylabel('10-fold Misclassification rate');
        %     xlabel('Learning cycle');

        %     Classifier = fitcknn(x_train,y_train,'NumNeighbors',5,'Standardize',1);
        %     Classifier = fitcecoc(x_train,y_train,'OptimizeHyperparameters','auto',...
        %     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
        %     'expected-improvement-plus'));
        %     Classifier = fitcecoc(x_train, y_train,'Coding','onevsone',...
        %         'Learners',svmtempate);
            tic;
            CPUSt = cputime();
            yfit = predict(Classifier,x_test);
            TestTime(Dir_i, Round, PCA_i) = toc;
            TestCPU(Dir_i, Round, PCA_i) = cputime() - CPUSt;
            
            RoundAccuracy(Round) = sum(yfit==y_test)/length(yfit);
            disp(num2str(RoundAccuracy(Round)));
            for i = 1 : length(y_test)
                confusion_matrices(y_test(i), yfit(i), Round) = ...
                    confusion_matrices(y_test(i), yfit(i), Round) + 1;
            end
        end
        disp(fullfile(FileDir, Dirs(Files(Dir_i)).name));
        disp(['PCA compoments =' num2str(num_for_PCA(PCA_i))]);
        All_Confusion{PCA_i} = confusion_matrices;
        All_Accuracy{PCA_i} = RoundAccuracy;
        All_Template{PCA_i} = ClassifierTemplate;
        Result = sum(confusion_matrices,3)
        RoundAccuracy
%         save(fullfile(Dir(end).folder,['Result-BandPass-AllPCA' SaveTime '.mat']),...
%         'All_Accuracy','All_Confusion','All_Template');
        save(fullfile(Dir(end).folder,['Result-Temp-AllPCA' SaveTime '.mat']),...
        'All_Accuracy','All_Confusion','All_Template');
        
    end
    catch ME
        if length(Files) == 1
            rethrow(ME);
        end
        Messages{ME_I,1} = getReport(ME);
        Messages_Idx{ME_I,1} = fullfile(FileDir, Dirs(Files(Dir_i)).name);
        if length(Files) == 1
            rethrow(ME);
        end
        fprintf('\n\n\n___________ERROR_START__________\n');
        disp(Messages_Idx{ME_I});
        disp(Messages{ME_I});
        fprintf('___________ERROR_END____________\n\n\n');
        ME_I = ME_I + 1;
    end
end
% tic
% tree = templateTree('Surrogate','On','MaxNumSplits',20);
% tree2 = templateTree('MinLeafSize',16);
% tree3 = templateTree('MinLeafSize',10);
% tree4 = templateTree('MinLeafSize',18);
% svmtempate = templateSVM('BoxConstraint', 38.666, 'KernelScale', 910.31);

% Result = sum(confusion_matrices,3)%;
% Accuracy = [diag(Result) ./ sum(Result,2);...
%     sum(diag(Result))/sum(sum(Result,2))]

%      [featrue_train,pstrain] = mapminmax(x_train');
%     % 将映射函数的范围参数分别置为0和4
%     pstrain.ymin = 0;
%     pstrain.ymax = 1;
%     % 对训练集进行[0,1]归一化
%    [featrue_train,pstrain] = mapminmax(featrue_train,pstrain);
%   % 测试数据处理
%     [feature_test,pstest] = mapminmax(x_test');
%     % 将映射函数的范围参数分别置为0和1
%     pstest.ymin = 0;
%     pstest.ymax = 1;
%     % 对测试集进行[0,1]归一化
%     [feature_test,pstest] = mapminmax(feature_test,pstest);
%     % 对训练集和测试集进行转置,以符合libsvm工具箱的数据格式要求
%     featrue_train = featrue_train';
%     feature_test = feature_test';
