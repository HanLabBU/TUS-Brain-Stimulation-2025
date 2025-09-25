%% transient and sustained/delayed vm
FS = 1000;

good_neurons = find(T_params.Good_ROI > 0 & T_params.Bad_Motion < 1);

mod_Vm_track_250ms = zeros(size(good_neurons));
mod_Vm_track_1000ms = zeros(size(good_neurons));

for count = 1:length(good_neurons)
    curr_neuron = good_neurons(count);

    preIdx = 1 * FS; 
    preWind = (preIdx - round(1 * FS) + 1):preIdx; 
    
    USWind_250ms = (preIdx):(preIdx + round(0.25 * FS)); % First 250 ms
    USWind_250to1000ms = (preIdx + round(0.25 * FS) + 1):(preIdx + round(1 * FS)); % 250–1000 ms

    Vm_mean_pre = squeeze(mean(vm_all(:, preWind, curr_neuron), 2, 'omitnan')); % Baseline
    Vm_mean_250ms = squeeze(mean(vm_all(:, USWind_250ms, curr_neuron), 2, 'omitnan')); % First 250 ms

    % --- Transient Period (0–250 ms) ---
    [p_250ms, ~, ~] = signrank(Vm_mean_250ms);
    if p_250ms < 0.05 && (mean(Vm_mean_250ms - Vm_mean_pre, 'omitnan')) > 0
        mod_Vm_track_250ms(count) = 1; % Significant modulation (increased)
    elseif p_250ms < 0.05 && (mean(Vm_mean_250ms - Vm_mean_pre, 'omitnan')) < 0
        mod_Vm_track_250ms(count) = -1; % Significant modulation (decreased)
    else
        mod_Vm_track_250ms(count) = 0; % No significant change
    end

    % --- Delayed Period (250–1000 ms) ---
    bin_size = round(0.25 * FS); % 250 ms bins
    num_bins = 3; % Three 250 ms bins within 250–1000 ms
    sustained_significant_bins = 0;

    for bin = 1:num_bins
        bin_start = preIdx + round(0.25 * FS) + (bin - 1) * bin_size + 1;
        bin_end = bin_start + bin_size - 1;

        Vm_mean_bin = squeeze(mean(vm_all(:, bin_start:bin_end, curr_neuron), 2, 'omitnan'));

        [p_bin, ~, ~] = signrank(Vm_mean_bin);
        if p_bin < 0.05 && mean(Vm_mean_bin - Vm_mean_pre, 'omitnan') > 0
            sustained_significant_bins = sustained_significant_bins + 1;
        end
    end

    if sustained_significant_bins >= 2
        % Check for two consecutive bins
        consecutive_bins = false;
        for bin = 1:(num_bins - 1)
            bin_start1 = preIdx + round(0.25 * FS) + (bin - 1) * bin_size + 1;
            bin_end1 = bin_start1 + bin_size - 1;

            bin_start2 = preIdx + round(0.25 * FS) + bin * bin_size + 1;
            bin_end2 = bin_start2 + bin_size - 1;

            Vm_mean_bin1 = squeeze(mean(vm_all(:, bin_start1:bin_end1, curr_neuron), 2, 'omitnan'));
            Vm_mean_bin2 = squeeze(mean(vm_all(:, bin_start2:bin_end2, curr_neuron), 2, 'omitnan'));

            [p_bin1, ~, ~] = signrank(Vm_mean_bin1);
            [p_bin2, ~, ~] = signrank(Vm_mean_bin2);

            if (p_bin1 < 0.05 && mean(Vm_mean_bin1 - Vm_mean_pre, 'omitnan') > 0) && ...
               (p_bin2 < 0.05 && mean(Vm_mean_bin2 - Vm_mean_pre, 'omitnan') > 0)
                consecutive_bins = true;
                break;
            end
        end

        if consecutive_bins
            if mod_Vm_track_250ms(count) ~= 0
                continue;
            else
                mod_Vm_track_1000ms(count) = 1; 
            end
        else
            mod_Vm_track_1000ms(count) = 0;
        end
    else
        mod_Vm_track_1000ms(count) = 0; 
    end
end

fprintf('\n%i/%i neurons 250ms modulated (increased)', sum(mod_Vm_track_250ms > 0), length(mod_Vm_track_250ms));
fprintf('\n%i/%i neurons 250ms modulated (decreased)', sum(mod_Vm_track_250ms < 0), length(mod_Vm_track_250ms));
fprintf('\n%i/%i neurons 250-1000ms modulated (increased)', sum(mod_Vm_track_1000ms > 0), length(mod_Vm_track_1000ms));
fprintf('\n%i/%i neurons 250-1000ms modulated (decreased)', sum(mod_Vm_track_1000ms < 0), length(mod_Vm_track_1000ms));

%% transient and sustained spiking
mod_Sp_track_250ms = zeros(size(good_neurons));
mod_Sp_track_1000ms = zeros(size(good_neurons));

for count = 1:length(good_neurons)
    curr_neuron = good_neurons(count);

    preIdx = 1 * FS;
    preWind = (preIdx - round(1 * FS) + 1):preIdx;

    USWind_250ms = (preIdx):(preIdx + round(0.25 * FS)); % First 250 ms
    USWind_250to1000ms = (preIdx + round(0.25 * FS) + 1):(preIdx + round(1 * FS)); % 250–1000 ms

    Vm_mean_pre = squeeze(mean(full_block_psth(:, preWind, curr_neuron), 2, 'omitnan')); % Baseline
    Vm_mean_250ms = squeeze(mean(full_block_psth(:, USWind_250ms, curr_neuron), 2, 'omitnan')); % First 250 ms

    % --- Transient Period (0–250 ms) ---
    [p_250ms, ~, ~] = signrank(Vm_mean_250ms);
    if p_250ms < 0.05 && (mean(Vm_mean_250ms - Vm_mean_pre, 'omitnan')) > 0
        mod_Sp_track_250ms(count) = 1;
    elseif p_250ms < 0.05 && (mean(Vm_mean_250ms - Vm_mean_pre, 'omitnan')) < 0
        mod_Sp_track_250ms(count) = -1; 
    else
        mod_Sp_track_250ms(count) = 0; 
    end

    % --- Sustained Period (250–1000 ms) ---
    bin_size = round(0.25 * FS); % 250 ms bins
    num_bins = 3; % Three 250 ms bins within 250–1000 ms
    sustained_significant_bins = 0;

    for bin = 1:num_bins
        bin_start = preIdx + round(0.25 * FS) + (bin - 1) * bin_size + 1;
        bin_end = bin_start + bin_size - 1;

        Vm_mean_bin = squeeze(mean(vm_all(:, bin_start:bin_end, curr_neuron), 2, 'omitnan'));

        [p_bin, ~, ~] = signrank(Vm_mean_bin);
        if p_bin < 0.05 && mean(Vm_mean_bin - Vm_mean_pre, 'omitnan') > 0
            sustained_significant_bins = sustained_significant_bins + 1;
        end
    end

    if sustained_significant_bins >= 2
        % Check for two consecutive bins
        consecutive_bins = false;
        for bin = 1:(num_bins - 1)
            bin_start1 = preIdx + round(0.25 * FS) + (bin - 1) * bin_size + 1;
            bin_end1 = bin_start1 + bin_size - 1;

            bin_start2 = preIdx + round(0.25 * FS) + bin * bin_size + 1;
            bin_end2 = bin_start2 + bin_size - 1;

            Vm_mean_bin1 = squeeze(mean(vm_all(:, bin_start1:bin_end1, curr_neuron), 2, 'omitnan'));
            Vm_mean_bin2 = squeeze(mean(vm_all(:, bin_start2:bin_end2, curr_neuron), 2, 'omitnan'));

            [p_bin1, ~, ~] = signrank(Vm_mean_bin1);
            [p_bin2, ~, ~] = signrank(Vm_mean_bin2);

            if (p_bin1 < 0.05 && mean(Vm_mean_bin1 - Vm_mean_pre, 'omitnan') > 0) && ...
               (p_bin2 < 0.05 && mean(Vm_mean_bin2 - Vm_mean_pre, 'omitnan') > 0)
                consecutive_bins = true;
                break;
            end
        end

        if consecutive_bins
            if mod_Sp_track_250ms(count) ~= 0
                continue;
            else
                mod_Sp_track_1000ms(count) = 1;
            end
        else
            mod_Sp_track_1000ms(count) = 0; 
        end
    else
        mod_Sp_track_1000ms(count) = 0; 
    end
end

fprintf('\n%i/%i neurons 250ms modulated (increased)', sum(mod_Sp_track_250ms > 0), length(mod_Sp_track_250ms));
fprintf('\n%i/%i neurons 250ms modulated (decreased)', sum(mod_Sp_track_250ms < 0), length(mod_Sp_track_250ms));
fprintf('\n%i/%i neurons 250-1000ms modulated (increased)', sum(mod_Sp_track_1000ms > 0), length(mod_Sp_track_1000ms));
fprintf('\n%i/%i neurons 250-1000ms modulated (decreased)', sum(mod_Sp_track_1000ms < 0), length(mod_Sp_track_1000ms));

%% adding to metadata

T_params.vm1000ms = mod_Vm_track_1000ms;
T_params.vm250ms = mod_Vm_track_250ms;

T_params.sp250ms = mod_Sp_track_250ms;
T_params.sp1000ms = mod_Sp_track_1000ms;

