% Train and evaluate the facial expression recognition
% clean the matlab
close all
clear
clc
% Options
file_path = fullfile('..', 'Data', 'Features.mat');
result_path = fullfile('..', 'Data', 'Result-%s-%s.mat');
kFold = 10;
users = [1]; % we can choose more than one user, it will process for each
process = {
    % BalanceDataset: balance the samples for all facial expressions, based
    % on the smallest number of samples: [50 44 48 50 50 50]->[44 44 44 44 44 44]
    'BalanceDataset', 0
    % Once: this is used to test to check whether the objective is good
    'Once', 0
    % RandomDataPrem: is used to shuffle the data
    'RandomDataPrem', 0
    % MapMinMax: map the variables into (0,1) for each feature
    'MapMinMax', 1
};
classifier_name = 'Ensemble'; % 'Ensemble', 'kNN', 'SVM'
% hyper-parameters for optimization
hyper_opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
        'MaxObjectiveEvaluations',30, 'UseParallel', false, ...
        'Kfold',kFold,'Verbose',1,'ShowPlots',false); 
% initialize
random_seed = 404;
process(:,2) = cellfun(@(x) {x==1},process(:,2));
opts = cell2struct(process(:,2),process(:,1),1);
load(file_path, 'FeatureData', 'Description', 'Expressions', 'FeatureName');
exp_num = length(Expressions);
user_num = length(users);
feature_num = size(FeatureData{1,1}, 1);
all_accuracy = zeros(user_num, kFold);
all_confumat = zeros(exp_num, exp_num, user_num);
all_template = cell(user_num);
save_time = datestr(datetime(),'mmddHHMMSS');
if opts.MapMinMax % map to (0, 1), for each user and each feature
    for ui = 1 : size(FeatureData,1)
        extremums = ones(feature_num, 2) .* [Inf, -Inf];
        for j = 1 : size(FeatureData,2)
            extremums(:,1) = min(extremums(:,1), min(FeatureData{ui,j},[],2));
            extremums(:,2) = max(extremums(:,2), max(FeatureData{ui,j},[],2));
        end
        for j = 1 : size(FeatureData,2)
            FeatureData{ui,j} = (FeatureData{ui,j} - extremums(:,1)) ./ ...
                (extremums(:,2) - extremums(:,1));
            FeatureData{ui,j}(isnan(FeatureData{ui,j})) = 0;
        end
    end
end
% for each user, we train a unique classifier to detect facial expressions
for ui = 1 : user_num
    rng(random_seed);
    u = users(ui);
    train_data = FeatureData(u,:)'; % just chose one user's data
    ActNum = cellfun(@(x) size(x,2), train_data); % count the sample number
    if opts.BalanceDataset  % make all expressions the same number
        ActNum(ActNum == 0) = [];
        Min = min(ActNum);
        for i = 1 : length(train_data)
            if ~isempty(train_data{i}) && size(train_data{i},2) > Min
                train_data{i}(:,Min+1:end) = [];
            end
        end
    end
    % we devide samples into kFold for each facial expression
    fold_group = zeros(kFold, exp_num); 
    train_label = cell(exp_num,1);
    for t = 1 : exp_num
        train_label{t} = ones(1,size(train_data{t},2));
        remainder = rem(length(train_label{t}), kFold);
        quotient = (length(train_label{t}) - remainder) / kFold;
        fold_group(:, t) = ones(kFold, 1) * quotient;
        fold_group(1:remainder, t) = fold_group(1:remainder, t) + 1;
        if opts.RandomDataPrem
            perm = randperm(size(train_data{t},2));
            train_data{t} = train_data{t}(:,perm);
        end
    end
    if opts.Once
        
        kFold = 1;
    end
    start_indices = ones(1, exp_num);
    for ri = 1 : kFold
        % Prepare the train and test sample:
        x_test = [];
        y_test = [];
        x_train = [];
        y_train = [];
        for ei = 1 : exp_num
            samples = train_data{ei}';
            s = start_indices(ei);
            t = s + fold_group(ri, ei) - 1;
            x_test = cat(1, x_test, samples(s:t,:));
            y_test = cat(1, y_test, ei * ones(t - s + 1,1));
            samples(s:t,:) = [];
            x_train= cat(1, x_train,samples);
            y_train = cat(1, y_train, ei * ones(size(samples, 1),1));
            start_indices(ei) = t + 1;
        end
        % choose classifier
        switch classifier_name
            case 'Ensemble'
                if ri == 1
                    template = fitcensemble(x_train,y_train,...
                        'OptimizeHyperparameters','auto',...
                        'HyperparameterOptimizationOptions', hyper_opts);
                    learner = template.ModelParameters.LearnerTemplates{1};
                end
                if strcmp(template.Method, 'Bag')
                    classifier = fitcensemble(x_train, y_train, 'Method',...
                        template.Method, 'Learners',learner,...
                        'NumLearningCycles',template.NumTrained);
                else
                    classifier = fitcensemble(x_train, y_train, 'Method',...
                        template.Method, 'Learners',learner,...
                        'NumLearningCycles',template.NumTrained,...
                        'LearnRate', template.ModelParameters.LearnRate);
                end
            case 'kNN'
                if ri == 1
                    template = fitcknn(x_train,y_train,...
                        'OptimizeHyperparameters','all',...
                        'HyperparameterOptimizationOptions', hyper_opts);
                    learner = template.ModelParameters;
                end
                classifier = fitcknn(x_train, y_train, ...
                    'NumNeighbors',template.NumNeighbors,...
                    'Distance',template.Distance,...
                    'DistanceWeight',template.DistanceWeight,...
                    'BreakTies', template.BreakTies,...
                    'NSMethod', template.NSMethod,...
                    'Exponent', learner.Exponent,...
                    'StandardizeData', learner.StandardizeData);
            case 'SVM'
                if ri == 1
                    template = fitcecoc(x_train,y_train,...
                        'OptimizeHyperparameters','all',...
                        'HyperparameterOptimizationOptions', hyper_opts);
                    learner = template.ModelParameters.BinaryLearners;
                end
                classifier = fitcecoc(x_train, y_train, 'Coding',...
                    template.CodingName, 'Learners',learner);
            otherwise
                error('Wrong classifier name');
        end
        % predict the test data that are not in the training set
        yfit = predict(classifier,x_test);
        all_accuracy(ui, ri) = sum(yfit==y_test)/length(yfit);
        disp(num2str(all_accuracy(ui, ri)));
        for i = 1 : length(y_test)
            all_confumat(y_test(i), yfit(i), ui) = ...
                all_confumat(y_test(i), yfit(i), ui) + 1;
        end
    end
    % save the result after each user are evaluated
        save(sprintf(result_path, classifier_name, save_time),...
            'all_accuracy', 'all_confumat');
end