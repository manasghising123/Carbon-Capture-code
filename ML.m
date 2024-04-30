clc


train_and_evaluate_svr('training_data.mat');

function train_and_evaluate_svr(filename)
    % Load training data
    data = load(filename);
    X_train = data.X_train;
    Y_train = data.Y_train;

    % Train SVR models for each output variable
    svr_models = cell(1, 5);
    for i = 1:5
        svr_models{i} = fitrsvm(X_train, Y_train(:, i), 'Standardize', true, 'KernelFunction', 'linear','OptimizeHyperparameters','all');
    end

    % Perform k-fold cross-validation for each model
    cv_models = cell(1, 5);
    for i = 1:5
        cv_models{i} = crossval(svr_models{i});
    end

    % Predictions from cross-validation for each model
    Y_cv_preds = cell(1, 5);
    for i = 1:5
        Y_cv_preds{i} = kfoldPredict(cv_models{i});
    end

    % Calculate correlation coefficients (R) between predicted and actual values for each model
    Rs = zeros(1, 5);
    for i = 1:5
        Rs(i) = corr(Y_train(:, i), Y_cv_preds{i});
    end

    % Plot actual vs predicted values for each output variable purity, recovery, productivity, energy_requirments
    labels = {'Purity', 'Recovery', 'Productivity', 'Energy Requirement', 'TCR'};
    
    
    for i = 1:5
        fig = figure;
        scatter(Y_train(:, i), Y_cv_preds{i});
        hold on;
        plot([min(Y_train(:, i)), max(Y_train(:, i))], [min(Y_train(:, i)), max(Y_train(:, i))], 'k--'); % Plot y = x line
        hold off;
        xlabel(['Actual ', labels{i}]);
        ylabel(['Predicted ', labels{i}]);
        title(['SVR Model Accuracy for ', labels{i}, ' (R = ', num2str(Rs(i)), ')']);
        grid on;
        saveas(fig, ['SVR_Model_Accuracy_', labels{i}, '.png']);
    end
    save('trained_svr_models.mat', 'svr_models')
end
