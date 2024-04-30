clc

% Generate training data if needed
% generate_training_data(500, 'training_data');

% Load training data
data = load('training_data.mat');
X_train = data.X_train;
Y_train = data.Y_train;

% Define neural network architecture
numInputs = size(X_train, 2); % Number of input features
numOutputs = size(Y_train, 2); % Number of output variables
hiddenLayerSize = 10; % Number of neurons in the hidden layer

% Create and train neural network models for each output variable
net_models = cell(1, numOutputs);
for i = 1:numOutputs
    net = fitnet(hiddenLayerSize);
    net = train(net, X_train', Y_train(:, i)');
    net_models{i} = net;
end

% Predictions from neural network models
Y_preds = cell(1, numOutputs);
for i = 1:numOutputs
    Y_preds{i} = net_models{i}(X_train');
    Y_preds{i} = Y_preds{i}';
end

% Calculate correlation coefficients (R) between predicted and actual values for each model
Rs = zeros(1, numOutputs);
for i = 1:numOutputs
    Rs(i) = corr(Y_train(:, i), Y_preds{i});
end

% Plot actual vs predicted values for each output variable
labels = {'Purity', 'Recovery', 'Productivity', 'Energy Requirement', 'TCR'};
for i = 1:numOutputs
    fig = figure;
    scatter(Y_train(:, i), Y_preds{i});
    hold on;
    plot([min(Y_train(:, i)), max(Y_train(:, i))], [min(Y_train(:, i)), max(Y_train(:, i))], 'k--'); % Plot y = x line
    hold off;
    xlabel(['Actual ', labels{i}]);
    ylabel(['Predicted ', labels{i}]);
    title(['Neural Network Accuracy for ', labels{i}, ' (R = ', num2str(Rs(i)), ')']);
    grid on;
    saveas(fig, ['Neural_Network_Model_Accuracy_', labels{i}, '.png']);
end

% Save trained neural network models
save('trained_nn_models.mat', 'net_models')
