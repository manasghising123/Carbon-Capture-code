clc;

%valid_samples_count = 0;

generate_training_data1(500, 'training_data');



function generate_training_data1(num_samples, filename)
    N = 10;
    
    SimParam = [
                    1169    -34500    -12000    28.5500
                    789     -24900    -5800     9.3000
                    782     -33000    -12000    9.3000
                    910     -31300    -18000    19.3000
                    427     -14000    -10000    2.8500
                    1200    -37853    -19434.6  17.5000
                    626     -25000    -12000    9.3000
                    1236    -32500    -18000    500
                    850     -32000    -12000    9.3000
                    1620    -55300    -12000    17.5000
                    1770    -42000    -19000    6.3300
                    1659    -33000    -12000    28.5500
                    950     -16100    -7000     29
                    1179    -27500    -10000    2.8000
                    1230    -30300    -12900    2.8000
                    1130    -36000    -15800    2
                ];

    IsothermPar = [
                        7.02257248438633	6.32616975723148	2.64746278252729e-12	2.26991391430361e-11	-33753.7010179947	-36418.7430295312	7.02257248438633	0	8.22811165712219e-09	0	-12000	0	0;
                        9.99992924661689	8.83878571667149	7.41699230433617e-11	1.29644046709677e-10	-24931.7081981987	-24931.7081981987	10	0	2.36861584011843e-08	0	-5800	0	0;
                        1	18.9200000000000	1.39584477212610e-09	5.91181315253409e-12	-33000	-33000	3.20000000000000	0	1.04013483089031e-08	0	-12000	0	0;
                        6.80000000000000	9.90000000000000	2.44000000000000e-11	1.39000000000000e-10	-42000	-24000	14	0	4.96000000000000e-10	0	-18000	0	0;
                        48	0	8.06000000000000e-10	0	-14000	-14000	48	0	1.11000000000000e-09	0	-10000	0	0;
                        6.21200000000000	7.15000000000000	4.37719502835270e-11	4.60879264890047e-13	-37853	-37853	11.9000000000000	0	1.38381081748143e-10	0	-19434.6000000000	0	0;
                        28.5538186436634	0	5.93852866269217e-11	0	-25000	-25000	28.5538186436634	0	1.22040239957074e-09	0	-12000	0	0;
                        1.53565167540578	2.67484603376329	5.85345193627719e-13	4.52325821519039e-12	-32500	-32500	1.53565167540578	0	5.01034516202895e-11	0	-18000	0	0;
                        6.81657823612286	0	8.43674415098552e-11	0	-32000	-32000	3.89059798361789	0	3.82759822483871e-09	0	-12000	0	0;
                        2.72727764341153	0.249599335971298	2.89127836156796e-12	2.00758727180358e-21	-52093.8779597488	-90236.5205072872	2.72727764341153	0	6.52000000000000e-09	0	-12000	0	0;
                        4.54093513692750	0	9.93128110818050e-13	0	-42000	-42000	2.18190720277024	0	1.53895154434375e-10	0	-19000	0	0;
                        5	3	9.46000000000000e-11	6.15000000000000e-16	-33000	-48000	12.7000000000000	0	4.29000000000000e-10	0	-12300	0	0;
                        7.26543818509888	0	1.21547436386518e-09	0	-16120.4419825087	-16120.4419825087	7.26543818509888	0	2.34589212168960e-10	0	-7020.02754740736	0	0;
                        7	5.10000000000000	3.20714028664748e-10	2.14894679105313e-11	-27500	-27500	6	0	5.18977441624449e-09	0	-10000	0	0;
                        6.57000000000000	3.13000000000000	1.44000000000000e-07	9.41000000000000e-07	-29243	-30800	9.60600000000000	0	7.35000000000000e-06	0	-12900	0	1;
                        3.09000000000000	2.54000000000000	8.65000000000000e-07	2.63000000000000e-08	-36641.2100000000	-35690.6600000000	5.84000000000000	0	2.50000000000000e-06	0	-15800	0	1
                    ];
    
    % Input parameter ranges (modify these according to your problem)
    L_range = [1.0, 7.0];           % Length of the column [m]
    P_0_range = [1.0e5, 10.0e5];    % Adsorption pressure [Pa]
    P_l_range = [0.1e5, 0.5e5];     % Evacuation pressure [Pa]
    t_ads_range = [10.0, 1000.0];   % Time of adsorption step [s]
    v_0_range = [0.1, 2];           % Inlet velocity of feed [m/s]
    alpha_range = [0.01, 0.99];     % Light product reflux ratio [-]
    beta_range = [0.0, 1.0];        % Heavy product reflux ratio [-]
    type_range = [1, 16];           % Type of adsorbent (integer values only)

    % Create Sobol' sequence object
    s = sobolset(8);  % 6-dimensional Sobol' sequence for 6 input parameters
    s = scramble(s,'MatousekAffineOwen');

    % Generate Sobol' sequence samples within the specified ranges
    samples = net(s, num_samples);

    % Scale the samples to the specified ranges
    scaled_samples = zeros(num_samples, 8);
    scaled_samples(:, 1) = samples(:, 1) * (L_range(2) - L_range(1)) + L_range(1);
    scaled_samples(:, 2) = samples(:, 2) * (P_0_range(2) - P_0_range(1)) + P_0_range(1);
    scaled_samples(:, 3) = samples(:, 3) * (P_l_range(2) - P_l_range(1)) + P_l_range(1);
    scaled_samples(:, 4) = samples(:, 4) * (t_ads_range(2) - t_ads_range(1)) + t_ads_range(1);
    scaled_samples(:, 5) = samples(:, 5) * (v_0_range(2) - v_0_range(1)) + v_0_range(1);
    scaled_samples(:, 6) = samples(:, 6) * (alpha_range(2) - alpha_range(1)) + alpha_range(1);
    scaled_samples(:, 7) = samples(:, 7) * (beta_range(2) - beta_range(1)) + beta_range(1);
    scaled_samples(:, 8) = samples(:, 8) * (type_range(2) - type_range(1)) + type_range(1);
    
    scaled_samples(:, 8) = round(scaled_samples(:, 8));
    
    scaled_samples(:, 8) = max(min(scaled_samples(:, 8), type_range(2)), type_range(1));

    % Initialize arrays to store input-output pairs
    X_train = zeros(num_samples, 8);                % Input parameters: L, P_0, P_l, t_ads, alpha, beta
    Y_train = zeros(num_samples, 5);  % Output variables: purity, recovery, productivity, energy_requirments, TCR

    % Generate training data using the scaled samples
    % Generate training data using the scaled samples
    for i = 1:num_samples
        try
            % Get the scaled sample
            scaled_sample = scaled_samples(i, :);

            % Extract input parameters
            L = scaled_sample(1);
            P_0 = scaled_sample(2);
            P_l = scaled_sample(3);
            t_ads = scaled_sample(4);
            v_0 = scaled_sample(5);
            alpha = scaled_sample(6);
            beta = scaled_sample(7);
            type = scaled_sample(8);

            % Call your PSA cycle simulation function to get the output variables
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, purity, recovery, productivity, energy_requirments, TCR] = PSACyclesample(L, P_0, P_l, t_ads, v_0, alpha, beta, N, type, SimParam, IsothermPar);

            % Store scalar output variables
            Y_train(i, :) = [purity, recovery, productivity, energy_requirments, TCR];

            % Store the input parameters only if the simulation is successful
            X_train(i, :) = [L, P_0, P_l, t_ads, v_0, alpha, beta, type];
            %valid_samples_count = valid_samples_count + 1
            disp(i)
        catch
            % If an error occurs, skip this sample and continue with the next
            % Skip adding this sample to X_train and Y_train
            X_train(i, :) = []; % Remove the corresponding row from X_train
            continue;
        end
    end

    % Remove any rows in X_train and Y_train that are all zeros (due to skipped samples)
    X_train(all(X_train == 0, 2), :) = [];
    Y_train(all(Y_train == 0, 2), :) = [];
    
    % Save the dataset as a MAT file
    save(filename + ".mat", 'X_train', 'Y_train');
end
