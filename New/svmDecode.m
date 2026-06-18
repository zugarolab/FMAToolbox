function [Accuracy, cvModel, testAccuracy, predictedFeatures, Model] = svmDecode(trainData, Features, varargin)

% svmDecode - Train SVM to decode features from data and evaluate accuracy
%
%  USAGE
%
%    [accuracy, model] = svmDecode(trainData, Features, <optional parameters>)
%
%  INPUTS
%     trainData     Training data matrix (n x m) : example : the average firing rate of m neurons for n trials
%     Features      Class labels / features to decode (n × 1) : example : the arm ID (1/2) for n trials
%
%    <optional parameters> can be :
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'KFold'        Number of cross-validation folds (default: 5) ; 5 means it will group features in 5 groups and will average  5 80% train/20% test combination
%     'testData'     Optional test data matrix (n2 x m), if not provided, accuracy will be from crossvalidated train data (if testFeatures not provided n2 must be equal to n1)
%     'testFeatures' Optional features to decode with test data (n2 × 1), if not provided, the accuracy will be computed by comparing with Features
%    =========================================================================
%
%  OUTPUTS
%     Accuracy          Cross-validated accuracy using training set
%     cvModel           Cross-validated SVM model   
%     testAccuracy      Accuracy using test set if provided (on train Features or on testFeatures if provided)
%     predictedFeatures Predicted Features with testData if provided
%     Model             SVM model (trained on full training set)
%
%  EXAMPLES
%
%     % Train model to decode arm identity based on neurons' firing rates
%     [acc, cvmodel, ~, ~, model] = svmDecode(FRs, armID) / predictedFeatures for each cross-validation set can be found in cvmodel.Trained{nCVset}.SupportVectorLabels  
%
%     % Train on mean FR of trial, try to predict with firing rate restricted to a part of the trial
%     [acc, cvmodel, acc2, predicted, model] = svmDecode(FRs, armID, 'testData', restrictedFR)
%
%     %  Train to decode arm identity based on neurons' firing rates on maze 1, try to predict arm identity in maze 2 based on neurons' firing rates on maze 2
%     [acc, cvmodel, acc2, predicted, model] = svmDecode(FRs, armID, 'testData', FRs2, 'testFeatures', armID2)

%% Optional parameters
p = inputParser;
addParameter(p, 'KFold', 2);
addParameter(p, 'testData', []);
addParameter(p, 'testFeatures', []);
parse(p, varargin{:});
KFold = p.Results.KFold;
testData = p.Results.testData;
testFeatures = p.Results.testFeatures;

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
        testAccuracy =  mean(predictedFeatures == Features);
    end
    
else 
    testAccuracy = [];
    predictedFeatures = [];
end


