clc

train_and_evaluate_et('training_data.mat');

function train_and_evaluate_et(filename)
    % Load training data
    data = load(filename);
    X_train = data.X_train;
    Y_train = data.Y_train;

    % Train Extra Trees models for each output variable
    et_models = cell(1, 5);
    for i = 1:5
        et_models{i} = TreeBagger(100, X_train, Y_train(:, i), 'Method', 'regression', 'MinLeafSize', 1, 'FBoot', 1);
    end

    % Perform k-fold cross-validation for each model
    k = 5; % Number of folds
    cv_indices = crossvalind('Kfold', size(X_train, 1), k); % Generate indices for k-fold cross-validation
    
    Y_cv_preds = cell(1, 5);
    for i = 1:k
        % Split data into training and validation sets
        train_indices = cv_indices ~= i;
        val_indices = cv_indices == i;
        X_train_fold = X_train(train_indices, :);
        Y_train_fold = Y_train(train_indices, :);
        X_val_fold = X_train(val_indices, :);
        
        % Train the model on training data
        et_model = TreeBagger(100, X_train_fold, Y_train_fold(:, i), 'Method', 'regression', 'MinLeafSize', 1, 'FBoot', 1);
        
        % Predict on validation data
        Y_cv_preds_fold = predict(et_model, X_val_fold);
        Y_cv_preds_fold = str2double(Y_cv_preds_fold); % Convert predictions to numeric
        
        % Store predictions for this fold
        Y_cv_preds{i} = Y_cv_preds_fold;
    end

    % Calculate correlation coefficients (R) between predicted and actual values for each model
    Rs = zeros(1, 5);
    for i = 1:5
        Rs(i) = corr(Y_train(:, i), cell2mat(Y_cv_preds(i)'));
    end

    % Plot actual vs predicted values for each output variable
    labels = {'Purity', 'Recovery', 'Productivity', 'Energy Requirement', 'TCR'};
    for i = 1:5
        fig = figure;
        scatter(Y_train(:, i), cell2mat(Y_cv_preds(i)'));
        hold on;
        plot([min(Y_train(:, i)), max(Y_train(:, i))], [min(Y_train(:, i)), max(Y_train(:, i))], 'k--'); % Plot y = x line
        hold off;
        xlabel(['Actual ', labels{i}]);
        ylabel(['Predicted ', labels{i}]);
        title(['Extra Trees Model Accuracy for ', labels{i}, ' (R = ', num2str(Rs(i)), ')']);
        grid on;
        saveas(fig, ['Extra_Trees_Model_Accuracy_', labels{i}, '.png']);
    end
    save('trained_et_models.mat', 'et_models')
end
