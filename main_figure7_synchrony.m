
%% Vm cross correlations at the neuron level -- showing Vm version, spiking is the same but with spiking data

FS = 1000;  % Sampling frequency in Hz
pre_period = 1:FS;  % Pre-stimulation period
post_period = (2*FS+1:3*FS);  % Post-stimulation period
time_window = 1000;  % Time window for synchrony (in ms)
total_time = length(pre_period);  % Total trial length in ms
max_lag = 500;  % Maximum lag (500 ms window for 1000-sample periods)

fov_list = T_meta.dir_name;  % Get FOV information from metadata

entrained_neurons = find([(T_params.vm250ms == 1) & (T_meta.Session ~= 81) & (T_meta.mouseID ~= 604104) & (T_meta.mouseID ~= 611293)]);

% Get unique FOVs containing entrained neurons
unique_fovs = unique(fov_list(entrained_neurons));
total_fovs = length(unique_fovs);

max_corr_pre = [];
max_corr_stim = [];
max_corr_post = [];

lag_at_max_corr_pre = [];
lag_at_max_corr_stim = [];
lag_at_max_corr_post = [];


% Loop through each FOV containing entrained neurons
for fov = 1:total_fovs
    neurons_in_fov = find(strcmp(fov_list, unique_fovs{fov}));
    total_neurons_in_fov = length(neurons_in_fov);
    
    if total_neurons_in_fov > 2
        fprintf('Processing FOV %d/%d (Entrained Neurons: %d)\n', fov, total_fovs, total_neurons_in_fov);
        lag_at_max_corr_pre_fov = [];
        lag_at_max_corr_stim_fov = [];
        lag_at_max_corr_post_fov = [];
        % Cross-correlation analysis for each pair of entrained neurons in the same FOV
        for i = 1:total_neurons_in_fov
            for j = i+1:total_neurons_in_fov
                neuron_i = neurons_in_fov(i);
                neuron_j = neurons_in_fov(j);

                num_trials = size(vm_all, 1);
                lag_at_max_corr_pre_trial = [];
                lag_at_max_corr_stim_trial = [];
                lag_at_max_corr_post_trial = [];
                for trial = 1:num_trials
                    % Extract voltage data for pre, stim, and post periods
                    vm_all = full_block_Vm;
                    V_data_i_pre = squeeze(vm_all(trial, pre_period, neuron_i));
                    V_data_j_pre = squeeze(vm_all(trial, pre_period, neuron_j));
    
                    V_data_i_stim = squeeze(vm_all(trial, pre_period+FS:pre_period+2*FS, neuron_i));
                    V_data_j_stim = squeeze(vm_all(trial, pre_period+FS:pre_period+2*FS, neuron_j));
    
                    V_data_i_post = squeeze(vm_all(trial, post_period, neuron_i));
                    V_data_j_post = squeeze(vm_all(trial, post_period, neuron_j));

                    V_data_i_pre(isnan(V_data_i_pre)) = 0;
                    V_data_j_pre(isnan(V_data_j_pre)) = 0;
                    V_data_i_stim(isnan(V_data_i_stim)) = 0;
                    V_data_j_stim(isnan(V_data_j_stim)) = 0;
                    V_data_i_post(isnan(V_data_i_post)) = 0;
                    V_data_j_post(isnan(V_data_j_post)) = 0;

                    if std(squeeze(vm_all(trial, :, neuron_i))) == 0 || std(squeeze(vm_all(trial, :, neuron_j))) == 0
                        continue;
                    end

                    [C_pre, lags_pre] = xcorr(V_data_i_pre, V_data_j_pre, max_lag, 'coeff');
                    [C_stim, lags_stim] = xcorr(V_data_i_stim, V_data_j_stim, max_lag, 'coeff');
                    [C_post, lags_post] = xcorr(V_data_i_post, V_data_j_post, max_lag, 'coeff');
    
                    % Find the index of zero lag
                    zero_lag_idx_pre = find(lags_pre == 0);
                    zero_lag_idx_stim = find(lags_stim == 0);
                    zero_lag_idx_post = find(lags_post == 0);
                    
                    % Extract and store the cross-correlation values at zero lag
                    corr_at_zero_lag_pre_trial = C_pre(zero_lag_idx_pre);
                    corr_at_zero_lag_stim_trial = C_stim(zero_lag_idx_stim);
                    corr_at_zero_lag_post_trial = C_post(zero_lag_idx_post);
                    
                    % Store the zero lag correlations in arrays
                    lag_at_max_corr_pre_trial = [lag_at_max_corr_pre_trial; corr_at_zero_lag_pre_trial];
                    lag_at_max_corr_stim_trial = [lag_at_max_corr_stim_trial; corr_at_zero_lag_stim_trial];
                    lag_at_max_corr_post_trial = [lag_at_max_corr_post_trial; corr_at_zero_lag_post_trial];
                end
                lag_at_max_corr_pre_fov = [lag_at_max_corr_pre_fov; nanmean(lag_at_max_corr_pre_trial)];
                lag_at_max_corr_stim_fov = [lag_at_max_corr_stim_fov; nanmean(lag_at_max_corr_stim_trial)];
                lag_at_max_corr_post_fov = [lag_at_max_corr_post_fov; nanmean(lag_at_max_corr_post_trial)];
            end
        end
        lag_at_max_corr_pre = [lag_at_max_corr_pre; nanmean(lag_at_max_corr_pre_fov)];
        lag_at_max_corr_stim = [lag_at_max_corr_stim; nanmean(lag_at_max_corr_stim_fov)];
        lag_at_max_corr_post = [lag_at_max_corr_post; nanmean(lag_at_max_corr_post_fov)];
    end
end

data_matrix = [lag_at_max_corr_pre, lag_at_max_corr_stim, lag_at_max_corr_post];
data_matrix_no_nan = data_matrix(~any(isnan(data_matrix), 2), :);  % Remove NaN rows

figure;
v = violinplot(data_matrix_no_nan, {'Pre-Stimulation', 'Stimulation', 'Post-Stimulation'});

% Turn off individual data points
for i = 1:length(v)
    v(i).ShowData = false;
end

ylabel('Cross-Correlation');
title('0 Lag Cross-Correlation');
ylim([0 1]); 

% Statistical Tests
alpha = 0.05;
num_comparisons = size(lag_at_max_corr_pre,1);  % Number of paired comparisons
alpha_corrected = alpha / num_comparisons;  % Bonferroni correction
data_matrix = [lag_at_max_corr_pre, lag_at_max_corr_stim, lag_at_max_corr_post];

[p_friedman, tbl, stats] = friedman(data_matrix, 1, 'off');
fprintf('Friedman Test p-value: %.5f\n', p_friedman);



%% signrank test

valid_rows = ~any(isnan(data_matrix), 2);  % remove any rows with NaN
data_clean = data_matrix(valid_rows, :);

% --- Wilcoxon signed-rank tests ---
[p_pre_stim, h_pre_stim] = signrank(data_clean(:,1), data_clean(:,2));
[p_stim_post, h_stim_post] = signrank(data_clean(:,2), data_clean(:,3));
[p_pre_post, h_pre_post] = signrank(data_clean(:,1), data_clean(:,3));

fprintf('Wilcoxon Signed-Rank Tests:\n');
fprintf('Pre vs Stim:  p = %.5f, h = %d\n', p_pre_stim, h_pre_stim);
fprintf('Stim vs Post: p = %.5f, h = %d\n', p_stim_post, h_stim_post);
fprintf('Pre vs Post:  p = %.5f, h = %d\n', p_pre_post, h_pre_post);