clc

%generate_training_data(500, 'training_data');

train_and_evaluate_rf('training_data.mat');

function train_and_evaluate_rf(filename)
    % Load training data
    data = load(filename);
    X_train = data.X_train;
    Y_train = data.Y_train;

    % Train Random Forest models for each output variable
    rf_models = cell(1, 5);
    for i = 1:5
        rf_models{i} = TreeBagger(50, X_train, Y_train(:, i), 'Method', 'regression');
    end

    % Predictions from cross-validation for each model
    Y_cv_preds = cell(1, 5);
    for i = 1:5
        Y_cv_preds{i} = predict(rf_models{i}, X_train);
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
        title(['Random Forest Model Accuracy for ', labels{i}, ' (R = ', num2str(Rs(i)), ')']);
        grid on;
        saveas(fig, ['Random_Forest_Model_Accuracy_', labels{i}, '.png']);
    end
    save('trained_rf_models.mat', 'rf_models')
    disp(tree_bagger_model);
end
