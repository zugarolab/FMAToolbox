function [predictedFeatures, accuracy, model] = svmDecode(trainingData, trainingFeatures, varargin)

% svmDecode - Train an SVM classifier to predict class labels / features from data and evaluate decoding accuracy.
%
%  USAGE
%
%    [predictedFeatures, accuracy, model] = svmDecode(trainingData, trainingFeatures, <optional parameters>)
%
%  INPUTS
%     trainingData      Training data matrix (n x m). Example: average firing rates of m neurons over n trials
%     trainingFeatures  Class labels / features to decode (n × 1) : example : the arm ID (1/2) for n trials
%
%    <optional parameters> can be :
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'testData'     Optional test data matrix (n2 x m) to predict features
%     'testFeatures' Optional features associated with the test data (n2 ×
%                    1 vector). If provided, accuracy will be computed based on how well
%                    the test features are predicted.
%     'kFold'        Number of cross-validation folds for the training set.
%                    Default: kFold=1 if testFeatures are provided (no cross-validation 
%                    needed for the training set, because the test set is separately 
%                    provided), or, if testFeatures are not provided,
%                    kFold=5 (accuracy evaluated using trainingFeatures).
%                    For example, with kFold = 5, the training data are
%                    split into five folds. The model is trained five times, 
%                    each time using 80% of the data for training and the
%                    remaining 20% for testing.
%    =========================================================================
%
%  OUTPUTS
%     predictedFeatures Predicted features with testData if provided
%     accuracy          Cross-validated (when kFold>1) accuracy using training set
%     model             SVM model (non-cross validated)
%
%  EXAMPLES
%
%     % Train SVM model to decode arm identity based on neurons' firing rates
%     [predictedFeatures, model, accuracy, ~, fullModel] = svmDecode(FRs, armID)
%
%     % Train on mean FR of trial, try to predict with firing rate restricted to a part of the trial
%     [acc, model, acc2, predicted, fullModel] = svmDecode(FRs, armID, 'testData', restrictedFR)
%
%     %  Train to decode arm identity based on neurons' firing rates on maze 1, try to predict arm identity in maze 2 based on neurons' firing rates on maze 2
%     [acc, model, acc2, predicted, fullModel] = svmDecode(FRs, armID, 'testData', FRs2, 'testFeatures', armID2)
%
% Copyright (C) 2026 by Théo Mathevet & Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%% Optional parameters

p = inputParser;
addParameter(p, 'kFold', []);
addParameter(p, 'testData', []);
addParameter(p, 'testFeatures', []);
parse(p, varargin{:});
kFold = p.Results.kFold;
testData = p.Results.testData;
testFeatures = p.Results.testFeatures;

if isempty(kFold), if isempty(testFeatures), kFold = 5; else, kFold = 1; end; end

% Train Support Vector Machine full model to decode features based on train Data
model = fitcsvm(trainingData, trainingFeatures, 'KernelFunction', 'linear', 'Standardize', true);

if kFold>1
    % Create kFold set for cross-validation
    cv = cvpartition(trainingFeatures, 'kFold', kFold);
    % Cross-validate model prediction
    cvmodel = crossval(model, 'CVPartition', cv);
else
    cvmodel = model;
end

if ~isempty(testData)
    predictedFeatures = predict(model, testData);
else
    if kFold>1
        predictedFeatures = kfoldPredict(cvmodel);
    else
        predictedFeatures = predict(model, trainingData);
    end
end

if ~isempty(testFeatures)
    accuracy = mean(predictedFeatures == testFeatures);
else % estimate accuracy using training data
    if kFold>1
        accuracy = 1 - kfoldLoss(cvmodel);
    else
        predictedTrainingFeatures = predict(model, trainingData);
        accuracy = mean(predictedTrainingFeatures == trainingFeatures);
    end
end


