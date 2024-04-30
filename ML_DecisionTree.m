clc

%generate_training_data(500, 'training_data');

train_and_evaluate_dt('training_data.mat');

function train_and_evaluate_dt(filename)
    % Load training data
    data = load(filename);
    X_train = data.X_train;
    Y_train = data.Y_train;

    % Train Decision Tree models for each output variable
    dt_models = cell(1, 5);
    for i = 1:5
        dt_models{i} = fitrtree(X_train, Y_train(:, i));
    end

    % Predictions from cross-validation for each model
    Y_cv_preds = cell(1, 5);
    for i = 1:5
        Y_cv_preds{i} = predict(dt_models{i}, X_train);
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
        title(['Decision Tree Accuracy for ', labels{i}, ' (R = ', num2str(Rs(i)), ')']);
        grid on;
        saveas(fig, ['Decision_Tree_Model_Accuracy_', labels{i}, '.png']);
    end
    save('trained_dt_models.mat', 'dt_models')
end
