function [Accuracy, cvModel, testAccuracy, Model] = svmDecode(trainData, Features, varargin)

% svmDecode - Train SVM to decode features from data and evaluate accuracy
%
%  USAGE
%
%    [accuracy, model] = svmDecode(trainData, Features, <optional parameters>)
%
%  INPUT
%     trainData     Training data matrix (n samples × n features)
%     Features      Class labels / features to decode (n samples × 1)
%
%    <optional parameters> can be :
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'KFold'        Number of cross-validation folds (default: 2)
%     'testData'     Optional test data matrix (n samples × n features), if not provided, accuracy will be from crossvalidated train data
%     'testFeatures' Optional features to decode with test data (n samples × 1), if not provided, the accuracy will be computed by comparing with Features
%    =========================================================================
%
%  OUTPUT
%     Accuracy      Cross-validated accuracy using training set
%     cvModel       Cross-validated SVM model
%     testAccuracy  Accuracy using test set (on train Features or on testFeatures if provided)
%     Model         SVM model (trained on full training set)
%
%  EXAMPLES
%
%     % Train model to decode arm identity based on neurons' firing rates
%     [acc, model] = svmDecode(FRs, armID)
%
%     % Train on mean FR of trial, try to predict with firing rate restricted to a part of the trial
%     [acc, model] = svmDecode(FRs, armID, 'testData', restrictedFR)
%
%     %  Train to decode arm identity based on neurons' firing rates on maze 1, try to predict arm identity in maze 2 based on neurons' firing rates on maze 2
%     [acc, model] = svmDecode(FRs, armID, 'testData', FRs2, 'testFeatures', armID2)

%% Optional parameters
p = inputParser;
addParameter(p, 'KFold', 2);
addParameter(p, 'testData', []);
addParameter(p, 'testFeatures', []);
parse(p, varargin{:});
KFold = p.Results.KFold;
testData = p.Results.testData;

%% Function

% Train Support Vector Machine model to decode features based on train Data
Model = fitcsvm(trainData, Features, 'KernelFunction', 'linear', 'Standardize', true);

% Create KFold set for cross-validation
cv = cvpartition(Features, 'KFold', KFold);

% Cross-validate model prediction
cvModel = crossval(Model, 'CVPartition', cv);

% Define cross-validated accuracy on train set
Accuracy = 1 - kfoldLoss(cvModel);

% If test set provided :
if ~isempty(testData)
    
    % Use model to predict features with test data 
    predictedFeatures = predict(Model, testData);
    
    % If test features are provided, compare to test features
    if ~isempty(testFeatures)
        testAccuracy = mean(predictedFeatures == testFeatures);
        
    % Else compare to train features    
    else
        testAccuracy = predictedFeatures;
    end

end


