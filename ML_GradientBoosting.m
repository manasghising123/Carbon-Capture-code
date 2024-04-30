clc

train_and_evaluate_gb('training_data.mat');

function train_and_evaluate_gb(filename)
    % Load training data
    data = load(filename);
    X_train = data.X_train;
    Y_train = data.Y_train;

    % Train Gradient Boosting models for each output variable
    gb_models = cell(1, 5);
    for i = 1:5
        gb_models{i} = fitensemble(X_train, Y_train(:, i), 'LSBoost', 100, 'Tree', 'Type', 'regression');
    end

    % Perform k-fold cross-validation for each model
    cv_models = cell(1, 5);
    for i = 1:5
        cv_models{i} = crossval(gb_models{i}, 'KFold', 5);
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

    % Plot actual vs predicted values for each output variable
    labels = {'Purity', 'Recovery', 'Productivity', 'Energy Requirement', 'TCR'};
    for i = 1:5
        fig = figure;
        scatter(Y_train(:, i), Y_cv_preds{i});
        hold on;
        plot([min(Y_train(:, i)), max(Y_train(:, i))], [min(Y_train(:, i)), max(Y_train(:, i))], 'k--'); % Plot y = x line
        hold off;
        xlabel(['Actual ', labels{i}]);
        ylabel(['Predicted ', labels{i}]);
        title(['Gradient Boosting Accuracy for ', labels{i}, ' (R = ', num2str(Rs(i)), ')']);
        grid on;
        saveas(fig, ['Gradient_Boosting_Model_Accuracy_', labels{i}, '.png']);
    end
    save('trained_gb_models.mat', 'gb_models')
end
