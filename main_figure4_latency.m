%% vm latency with baseline normalization 

FS = 1000;  % Sampling rate (Hz)
preIdx = 1 * FS;  % 1 second pre-stim baseline
stim_dur = 1 * FS; % 1 second stim
time_points = (-preIdx:stim_dur*2) / FS * 1000; % Time in ms

% transiently modulated neurons
transient_neurons = find( T_meta.PRF ~= 40 & ...
                         (T_params.vm250ms == 1) & ...
                         (T_meta.Session ~= 81) & ...
                         (T_meta.mouseID ~= 604104) & ...
                         (T_meta.mouseID ~= 611293));

num_neurons = length(transient_neurons);
p_values = NaN(num_neurons, stim_dur); 
first_significant_times = NaN(num_neurons, 1);

for n = 1:num_neurons
    curr_neuron = transient_neurons(n);
    Vm_trials = squeeze(full_block_trace(:, :, curr_neuron)); 
    num_trials = size(Vm_trials, 1);

    % Z-score each trial
    Vm_trials_z = NaN(size(Vm_trials));
    for trial = 1:num_trials
        base = Vm_trials(trial, 1:preIdx);
        mu = mean(base, 'omitnan');
        sigma = std(base, 'omitnan');
        if sigma == 0, sigma = 1; end
        Vm_trials_z(trial, :) = (Vm_trials(trial, :) - mu) / sigma;
    end

    % Get stim period (1s after preIdx)
    stim_Vm_z = Vm_trials_z(:, preIdx+1 : preIdx+stim_dur); 

    % Compute p-values vs. zero baseline
    for t = 1:stim_dur
        p_values(n, t) = ranksum(stim_Vm_z(:, t), zeros(num_trials, 1));
    end

    % Find first occurrence of consecutive significant p-values
    sig_mask = p_values(n, :) < 0.05; 
    sig_mask = [0, sig_mask]; % pad to catch start
    diff_mask = diff(sig_mask);
    run_starts = find(diff_mask == 1); % where a run starts
    run_lengths = diff(find([diff_mask 1])); % lengths of runs

    valid_runs = find(run_lengths >= 1);
    if ~isempty(valid_runs)
        first_sig_idx = run_starts(valid_runs(1)); % index into stim
        first_significant_times(n) = first_sig_idx; % in samples
    end
end

% Convert first significant sample indices to ms
valid_idx = ~isnan(first_significant_times);
first_significant_times_ms = (first_significant_times(valid_idx) / FS) * 1000;

% Stats
mean_latency = mean(first_significant_times_ms);
sem_latency = std(first_significant_times_ms) / sqrt(sum(valid_idx));

fprintf('\nMean first significant Vm latency (Z-scored) across neurons: %.2f ± %.2f ms (n = %d neurons, %d - %d ms, %d < 10 ms)\n', ...
    mean_latency, sem_latency, sum(valid_idx), min(first_significant_times), max(first_significant_times), sum(first_significant_times < 10));

figure;
histogram(first_significant_times_ms, 'BinWidth', 2, 'FaceColor', [0.2 0.4 0.6], 'EdgeColor', 'k');
xlabel('Latency to First Significant Vm Modulation (ms)');
ylabel('Number of Neurons');
title(sprintf('Distribution of First Significant Vm Responses (n = %d)', sum(valid_idx)));
xlim([0, 60]);
grid on;