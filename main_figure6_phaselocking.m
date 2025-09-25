%% Phase-Locking and Entrainment Analysis

FS = 1000; % sampling frequency
preIdx = 1*FS; 
postIdx = 1*FS;
stim_inds = (preIdx+1):(preIdx+postIdx);
stim_time = (1:length(stim_inds))/FS;

bl_inds = 1:preIdx;
bl_time = (1:length(bl_inds))/FS;

post_inds = (preIdx+postIdx+1):(size(full_block_Vm,2)-1);
post_time = (1:length(post_inds))/FS;

PRFs_loop = sort(unique(T_params.PRF),'ascend');

color_map = [0.8 0.3 0.3;...
            0.3 0.3 0.8];

% set filter/hilbert params
Fn = FS/2; % nyquist freq
frs = 1:50; % Hz
bandwidth_factor = 0.05; % Bandwidth (±5% around the PRF)

% spike shift for spike PLV
sp_shift = 0;

entrain_thresh_all = zeros(1,length(good_neurons));
entrain_mod_all = zeros(1,length(good_neurons));
entrain_shuff = cell(1,length(good_neurons));
entrain_PLV = zeros(1,length(good_neurons));
allFreqPLV = [];

phase_angles_all_BL = cell(length(good_neurons), 1); % Baseline
phase_angles_all_US = cell(length(good_neurons), 1); % Stimulation
phase_angles_all_Po = cell(length(good_neurons), 1); % Post-stimulation


for PRF_ind = 1:length(PRFs_loop)
    fprintf('\nBeginning PRF %i\n',PRFs_loop(PRF_ind))
    PRF = PRFs_loop(PRF_ind); % select PRF cells

    cells_curr_PRF = find(T_params.PRF == PRF);

    % Dynamically adjust filtering range based on PRF
    filter_low = PRF * (1 - bandwidth_factor);  
    filter_high = PRF * (1 + bandwidth_factor); 
    
    filter_low = max(1, filter_low);  
    filter_high = min(Fn, filter_high);  
    
    fprintf('Applying bandpass filter: %.2f - %.2f Hz around PRF: %.2f Hz\n', filter_low, filter_high, PRF);
    
    PLV_stim_PRF_BL = [];
    PLV_Vm_PRF_BL = [];
    PLV_Vm_Stim_PRF_BL = [];
    PLV_Vm_PRF_BL2 = [];
    
    PLV_stim_PRF_US = [];
    PLV_Vm_PRF_US = [];
    PLV_Vm_PRF_US2 = [];
    
    PLV_Vm_Stim_PRF_US = [];
    PLV_Vm_Stim_PRF_US2 = [];
    
    PLV_stim_PRF_Po = [];
    PLV_Vm_PRF_Po2 = [];
    PLV_Vm_Stim_PRF_Po = [];
    
    neuron_ids = [];
    for cell_ind = 1:length(cells_curr_PRF) 
        vm_all = full_block_trace;
        curr_cell = cells_curr_PRF(cell_ind);        
        num_trials_curr = sum(~isnan(vm_all(:, 500, curr_cell)));
        Vm_curr_cell = squeeze(vm_all(1:num_trials_curr, :, curr_cell)); 
        Sp_curr_cell = squeeze(full_block_raster(1:num_trials_curr, :, curr_cell)); 
        Sp_curr_cell(isnan(Sp_curr_cell)) = 0; 
        curr_pulses = pulses_save_all{curr_cell}; % US pulse idx for current cell
        
        pulses_curr_cell = zeros(size(Sp_curr_cell));
        pulses_curr_cell(:,curr_pulses) = 1; % make pseudo-raster for US pulse onset times
        
        data_all = [];
        wavD_Vm = zeros(size(Vm_curr_cell,2),length(frs),num_trials_curr);
        for trial_ind = 1:num_trials_curr % loop through trials
            Vmtr = Vm_curr_cell(trial_ind,:);
            Vmtr(isnan(Vmtr)) = 0;

            % filter trial Vm trace   
            filt_Vm = zeros(length(Vmtr),length(frs));
            for steps = 1:length(frs)
                FB = [frs(steps)*0.95 frs(steps)*1.05];                
                [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                filt_Vm(:,steps) = hilbert(filtfilt(B,A,Vmtr));
            end
            
            % Compile filtered signal wavelet across trials
            wavD_Vm(:,:,trial_ind) = filt_Vm;

        end

        wavD_Vm_post = zeros(1000,length(frs),num_trials_curr);
        for trial_ind = 1:num_trials_curr % loop through trials
            Vmtr = Vm_curr_cell(trial_ind,1001:2000);
            Vmtr(isnan(Vmtr)) = 0;

            % filter trial Vm trace
            filt_Vm = zeros(length(Vmtr),length(frs));
            for steps = 1:length(frs)
                FB = [frs(steps)*0.95 frs(steps)*1.05];                
                [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                filt_Vm(:,steps) = hilbert(filtfilt(B,A,Vmtr));
            end
            
            wavD_Vm_post(:,:,trial_ind) = filt_Vm;

        end

        all_Sp = Sp_curr_cell';
        all_St = pulses_curr_cell';
       
        % Get phase of wavelet 
        M_vm = angle(wavD_Vm);
        M_vm_post = angle(wavD_Vm_post);
        
        phase_angles_trials_BL = cell(1, num_trials_curr);
        phase_angles_trials_US = cell(1, num_trials_curr);
        phase_angles_trials_Po = cell(1, num_trials_curr);
        
        for trial_ind = 1:num_trials_curr
            stim_points_US = find(all_St(stim_inds, trial_ind));
            
            stim_points_BL = stim_points_US - preIdx; 
            stim_points_Po = stim_points_US + postIdx; 
            
            stim_points_BL(stim_points_BL <= 0) = [];
            stim_points_Po(stim_points_Po > length(post_inds)) = [];
            
            phase_angles_trials_BL{trial_ind} = angle(wavD_Vm(stim_inds-1000, PRF, trial_ind)); % Baseline
            phase_angles_trials_US{trial_ind} = angle(wavD_Vm(stim_inds, PRF, trial_ind)); % Stimulation
            phase_angles_trials_Po{trial_ind} = angle(wavD_Vm(stim_inds+1000, PRF, trial_ind)); % Post-stimulation
        end
        
        phase_angles_all_BL{curr_cell} = phase_angles_trials_BL;
        phase_angles_all_US{curr_cell} = phase_angles_trials_US;
        phase_angles_all_Po{curr_cell} = phase_angles_trials_Po;

        % Get Spike-Vm phase locking value across time inds
        [PLV_Vm_BL, PLV_Vm_BL2] = phase_locking_value(M_vm(bl_inds,:,:), all_St(stim_inds,:), sp_shift);              
        [PLV_Vm_US, PLV_Vm_US2] = phase_locking_value(M_vm(stim_inds,:,:), all_St(stim_inds,:), sp_shift);              
        [PLV_Vm_Po, PLV_Vm_Po2] = phase_locking_value(M_vm(post_inds,:,:), all_St(stim_inds,:), sp_shift);                       
        
        % Get Vm-Stim phase locking value    
        all_St2 = repmat(all_St(stim_inds,:), 2, 1);
        [PLV_StimVm_US, PLV_StimVm_US2] = phase_locking_value(M_vm(stim_inds,:,:), all_St(stim_inds,:), sp_shift);      

        % Initialize shuffled PLV array
        iter = 1500; % Number of shuffles
        PLV_Vm_Stim_PRF_sh = [];
        
        M_vm_selected = M_vm([bl_inds post_inds], :);

        % Shuffle and collect PLV values
        for count = 1:iter
            stim_shift_rand = randi([1, 1500]); 
            [PLV_StimVm_sh, PLV_StimVm_sh2] = phase_locking_value(M_vm([bl_inds stim_inds post_inds],frs == PRF,:),...
                all_St([bl_inds stim_inds post_inds],:), stim_shift_rand); 
            PLV_Vm_Stim_PRF_sh = [PLV_Vm_Stim_PRF_sh, PLV_StimVm_sh]; % Collect shuffled PLV values
        end

        entrain_thresh_curr = prctile(PLV_Vm_Stim_PRF_sh(1,:), 95); % 95th percentile of shuffled baseline PLV
        
        entrain_PLV(curr_cell) = PLV_StimVm_US(frs == PRF);
        entrain_thresh_all(curr_cell) = entrain_thresh_curr;

        if entrain_thresh_curr < entrain_PLV(curr_cell)
            entrain_mod_all(curr_cell) = 1; % significantly entrained phase to stim 
        end

        entrain_shuff{curr_cell} = PLV_Vm_Stim_PRF_sh;
     
        PLV_Vm_PRF_BL2 = [PLV_Vm_PRF_BL2 PLV_Vm_BL2'];
        PLV_Vm_PRF_US2 = [PLV_Vm_PRF_US2 PLV_Vm_US2'];
        PLV_Vm_PRF_Po2 = [PLV_Vm_PRF_Po2 PLV_Vm_Po2'];

        PLV_Vm_Stim_PRF_US2 = [PLV_Vm_Stim_PRF_US2 PLV_StimVm_US2];   
        PLV_Vm_Stim_PRF_US = [PLV_Vm_Stim_PRF_US PLV_StimVm_US];  

        T_meta.PLV(curr_cell) = entrain_PLV(curr_cell);
        T_meta.entrained(curr_cell) = entrain_mod_all(curr_cell);
        allFreqPLV = [allFreqPLV PLV_StimVm_US];
        fprintf('Neuron %i/%i complete\n',cell_ind,length(cells_curr_PRF))

    end
end

%% inter-trial coherence (all)

FS = 1000; % Sampling frequency (Hz)
Fn = FS / 2; % Nyquist frequency
frequencies = 1:50; % Frequency range to analyze for ITC (1-20 Hz)
time_window = 3; % Time window in seconds for analysis (adjust as needed)
num_trials = sum(~isnan(full_block_raster(:,500,:)), 1); % Get number of trials per neuron

window_size = 0.05; % Window size in seconds for wavelet analysis

itc_10hz = zeros(length(frequencies), FS * time_window + 1); % ITC for 10 Hz
itc_40hz = zeros(length(frequencies), FS * time_window + 1); % ITC for 40 Hz

neurons_10hz = find(T_meta.PRF == 10 & T_meta.entrained == 1); % 10 Hz neurons
neurons_40hz = find(T_meta.PRF == 40 & T_meta.entrained == 1); % 40 Hz neurons

itc_10hz = compute_itc(neurons_10hz, vm_all, frequencies, num_trials, time_window);
itc_40hz = compute_itc(neurons_40hz, vm_all, frequencies, num_trials, time_window);

time_vector = linspace(0, time_window, FS * time_window + 1);

Y = 1000 * [4.1962 3.5087 2.8783 2.3026 1.7775 1.2967 0.8516 0.4318]';
Y = interp1(Y, linspace(1, 8))';

Y_normalized = (Y - min(Y)) / (max(Y) - min(Y));
custom_cmap = flipud(interp1(linspace(0, 1, length(Y_normalized)), parula(length(Y_normalized)), Y_normalized));

figure(1);

% Plot for 10 Hz
subplot(2, 2, 1); % Main plot for 10 Hz
imagesc(time_vector, frequencies, itc_10hz); 
axis xy; 
colormap(parula); 
colorbar; 
xlabel('Time from Stim Onset (s)'); 
ylabel('Frequency (Hz)'); 
title('Inter-Trial Coherence (ITC) for 10 Hz Entrained Neurons'); 
ylabel(colorbar, 'ITC'); 
caxis([.2, .4]);

subplot(2, 2, 3); % Zoomed-in plot for 10 Hz
zoomed_freqs_10hz = frequencies(frequencies >= 5 & frequencies <= 15);
zoomed_times_10hz = time_vector(time_vector >= 0.5 & time_vector <= 2.5);
zoomed_itc_10hz = itc_10hz(frequencies >= 5 & frequencies <= 15, time_vector >= 0.5 & time_vector <= 2.5);
imagesc(zoomed_times_10hz, zoomed_freqs_10hz, zoomed_itc_10hz);
axis xy; 
colormap(parula); 
colorbar; 
xlabel('Time from Stim Onset (s)'); 
ylabel('Frequency (Hz)'); 
title('Zoomed ITC: 10 Hz (5-15 Hz, -0.5 to 1.5 s)');
caxis([.2, .4]);

% Plot for 40 Hz
subplot(2, 2, 2); % Main plot for 40 Hz
imagesc(time_vector, frequencies, itc_40hz); 
axis xy; 
colormap(parula); 
colorbar; 
xlabel('Time from Stim Onset (s)'); 
ylabel('Frequency (Hz)'); 
title('Inter-Trial Coherence (ITC) for 40 Hz Entrained Neurons'); 
ylabel(colorbar, 'ITC'); 
caxis([.1, .5]);

subplot(2, 2, 4); % Zoomed-in plot for 40 Hz
zoomed_freqs_40hz = frequencies(frequencies >= 35 & frequencies <= 45);
zoomed_times_40hz = time_vector(time_vector >= .5 & time_vector <= 2.5);
zoomed_itc_40hz = itc_40hz(frequencies >= 35 & frequencies <= 45, time_vector >= 0.5 & time_vector <= 2.5);
imagesc(zoomed_times_40hz, zoomed_freqs_40hz, zoomed_itc_40hz);
axis xy; 
colormap(parula); 
colorbar; 
xlabel('Time from Stim Onset (s)'); 
ylabel('Frequency (Hz)'); 
title('Zoomed ITC: 40 Hz (35-45 Hz, -0.5 to 1.5 s)');
caxis([.1, .5]);

%% ITC band plotted
FS = 1000; % Sampling frequency (Hz)
time_window = 3; % Time window in seconds for analysis
frequencies = 1:50; % Frequency range to analyze for ITC
time_vector = linspace(-1, 2, FS * time_window + 1); % Time vector for plotting

freq_idx_10hz = find(frequencies == 10); % Index for 10 Hz
freq_idx_40hz = find(frequencies == 40); % Index for 40 Hz

num_trials = sum(~isnan(full_block_raster(:, 500, :)), 1);

neurons_10hz = find(T_meta.PRF == 10 & T_meta.entrained > 0); % 10 Hz neurons
neurons_40hz = find(T_meta.PRF == 40 & T_meta.entrained > 0); % 40 Hz neurons

[mean_itc_10hz, sem_itc_10hz] = compute_itc_sem(neurons_10hz, vm_all, frequencies, num_trials, time_window);
[mean_itc_40hz, sem_itc_40hz] = compute_itc_sem(neurons_40hz, vm_all, frequencies, num_trials, time_window);

baseline_period = time_vector < 0;

% Compute the mean baseline ITC for 10 Hz and 40 Hz
baseline_mean_10hz = mean(mean_itc_10hz(10, baseline_period));
baseline_mean_40hz = mean(mean_itc_40hz(40, baseline_period));

figure;
subplot(1, 2, 1);
hold on;
plot(time_vector, mean_itc_10hz(10, :), 'b', 'LineWidth', 2); % Mean ITC for 10 Hz
fill([time_vector, fliplr(time_vector)], ...
     [mean_itc_10hz(10, :) + sem_itc_10hz(10, :), fliplr(mean_itc_10hz(10, :) - sem_itc_10hz(10, :))], ...
     'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % SEM shading
yline(baseline_mean_10hz, 'k--', 'LineWidth', 1.5); % Baseline mean as a dashed line
xlabel('Time (s)');
ylabel('ITC');
title('ITC for 10 Hz Neurons');
hold off;

% Plot for 40 Hz
subplot(1, 2, 2);
hold on;
plot(time_vector, mean_itc_40hz(40, :), 'r', 'LineWidth', 2); % Mean ITC for 40 Hz
fill([time_vector, fliplr(time_vector)], ...
     [mean_itc_40hz(40, :) + sem_itc_40hz(40, :), fliplr(mean_itc_40hz(40, :) - sem_itc_40hz(40, :))], ...
     'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % SEM shading
yline(baseline_mean_40hz, 'k--', 'LineWidth', 1.5); % Baseline mean as a dashed line
xlabel('Time (s)');
ylabel('ITC');
title('ITC for 40 Hz Neurons');
hold off;


%% average across PRFs: rose plots, Vm trace, heatmap
PRF_vals = [10, 40];
colors = {[147, 33, 31]/255, [5, 93, 113]/255}; % red = 10 Hz, blue = 40 Hz
FS = 1000;

for i = 1:length(PRF_vals)
    PRF = PRF_vals(i);
    curr_color = colors{i};
    entrained_cells = find(T_meta.entrained == 1 & T_params.PRF == PRF);

    neuron_avg_angles = [];
    Vm_z_all = [];  % stores normalized Vm per neuron

    for c = 1:length(entrained_cells)
        cell_idx = entrained_cells(c);
        trial_phases = phase_angles_all_US{cell_idx}; 
        pulse_trials = pulses_save_all{cell_idx};     

        angles_this_neuron = [];
        Vm_this_neuron = squeeze(vm_all(:, :, cell_idx));
        n_trials = sum(~isnan(Vm_this_neuron(:, 1000)));
        
        if n_trials == 0, continue; end
        
        Vm_mean = mean(Vm_this_neuron(1:n_trials,:), 1, 'omitnan');
        baseline = 1:FS;  % first second
        mu = mean(Vm_mean(baseline), 'omitnan');
        sigma = std(Vm_mean(baseline), [], 'omitnan');

        Vm_z = (Vm_mean - mu) / sigma;
        Vm_z_all = [Vm_z_all; Vm_z];

        for trial_num = 1:length(trial_phases)
            if isempty(trial_phases{trial_num}) || trial_num > length(pulse_trials)
                continue
            end
            curr_phases = trial_phases{trial_num}; 
            curr_pulses = pulse_trials - 1000;

            valid_inds = curr_pulses(curr_pulses > 0 & curr_pulses <= length(curr_phases));
            if ~isempty(valid_inds)
                angles_this_neuron = [angles_this_neuron; curr_phases(valid_inds)'];
            end
        end

        if ~isempty(angles_this_neuron)
            avg_angle = angle(mean(exp(1i * angles_this_neuron))); % same as before
            neuron_avg_angles = [neuron_avg_angles; avg_angle];
        end
    end

    % === Sort heatmap by mean Vm from 1001 to 1250 ms ===
    if ~isempty(Vm_z_all)
        sort_window = 1001:1250;
        Vm_sort_vals = mean(Vm_z_all(:, sort_window), 2, 'omitnan');
        [~, sort_idx] = sort(Vm_sort_vals, 'descend');
        Vm_z_all_sorted = Vm_z_all(sort_idx, :);
    else
        Vm_z_all_sorted = [];
    end

    % --- Plot Summary Figure ---
    figure('Name', sprintf('%d Hz Entrained Summary', PRF), 'Position', [100 100 1400 400]);
    tiledlayout(1, 3);

    % 1. Polar histogram of average phase angles
    nexttile(1);
    bin_edges = deg2rad(0:20:360);
    mean_angle = circ_mean(neuron_avg_angles(:));
    r = circ_r(neuron_avg_angles(:));
    [pval, ~] = circ_rtest(neuron_avg_angles(:));
    pval = min(pval * length(neuron_avg_angles(:)), 1) % Bonferroni

    polarhistogram(neuron_avg_angles, bin_edges, ...
        'FaceColor', curr_color, ...
        'EdgeColor', 'k', ...
        'FaceAlpha', 0.7, ...
        'Normalization', 'probability');
    title(sprintf('%d Hz Entrained\np = %.4f, r = %.2f', PRF, pval, r));

    % 2. Heatmap of z-scored Vm (sorted)
    nexttile(2);
    imagesc(Vm_z_all_sorted);
    colormap_matrix = flipud(cbrewer('div', 'RdBu', 500));
    colormap_matrix(colormap_matrix > 1) = 1;
    colormap_matrix(colormap_matrix < 0) = 0;
    colormap(colormap_matrix);
    caxis([-5, 5]);
    colorbar;
    xline(1000, '-k');
    xlabel('Time (ms)');
    ylabel('Neuron #');
    title('Z-scored Vm (Sorted by 1001–1250 ms)');

    % 3. Average z-scored Vm trace with SEM
    nexttile(3);
    mean_z = mean(Vm_z_all, 1, 'omitnan');
    sem_z = std(Vm_z_all, [], 1, 'omitnan') ./ sqrt(size(Vm_z_all, 1));
    x = (1:length(mean_z)) / FS - 1;  % centered on stim

    fill([x, fliplr(x)], ...
         [mean_z + sem_z, fliplr(mean_z - sem_z)], ...
         curr_color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hold on;
    plot(x, mean_z, 'Color', curr_color, 'LineWidth', 2);
    xline(0, '-k');
    xlabel('Time (s)');
    ylabel('Z-scored Vm');
    title('Avg Normalized Vm Across Entrained Neurons');
    axis tight;

    % --- Calculate descriptive stats
    std_angle = circ_std(neuron_avg_angles(:));
    n_neurons = length(neuron_avg_angles(:));
    
    fprintf('Across all %d Hz entrained neurons, the preferred phase angles formed a tight cluster around %.1f ± %.1f° (mean ± std, n=%d pulses).\n', ...
        PRF, rad2deg(mean_angle), rad2deg(std_angle), n_neurons);

end


%% functions

function [PLV_output, PLV_output2] = phase_locking_value(wavD, spikes, timshift)  
    [t, frs, tr] = size(wavD);   
    phase_angles = {};

    mat = [];    
    for trial_ind = 1:tr    
        M = exp(1i.*((squeeze(wavD(:,:,trial_ind))))); 
        s = find(spikes(:,trial_ind));
        s = s-timshift; 
        s(s<=0) = [];    
        mat = [mat, M(s,:)']; % The normal PLV formula  
        phase_angles{trial_ind} = angle(M(s,:)); 
    end    
        
    NT = sum(~isnan(mat(1,:)));    
    Z = abs(mean(mat,2,'omitnan'));  % MVL    
    T = Z.^2;    
    PLV_output = (((1/(NT-1))*((T.*NT-1))));  % adjusted MLV (PPC)    
%     PLV_output2 = mat; % normal PLV
    PLV_output(PLV_output < 0) = 0;
    PLV_output2 = Z;
    
    % % Exclusion criteria for how many events were detected
    % if NT <= 10    
    %     PLV_output = PLV_output.*NaN;    
    % end
end


function itc_values = compute_itc(neurons, vm_all, frequencies, num_trials, time_window)
    FS = 1000;
    Fn = FS / 2;
    itc_values = zeros(length(frequencies), FS * time_window + 1);
    
    for cell_ind = 1:length(neurons)
        curr_cell = neurons(cell_ind); % Current neuron index
        
        Vm_curr_cell = squeeze(vm_all(:, :, curr_cell));
        Vm_curr_cell(isnan(Vm_curr_cell)) = 0; % Replace NaNs with zero

        wavelet_matrix = zeros(length(frequencies), FS * time_window + 1, num_trials(cell_ind));

        for trial_ind = 1:num_trials(cell_ind)
            Vm_trace = Vm_curr_cell(trial_ind, :); % Extract voltage trace for this trial
            Vm_trace(isnan(Vm_trace)) = 0; % Replace NaNs with zero
            
            for freq_ind = 1:length(frequencies)
                fb = [frequencies(freq_ind) * 0.95, frequencies(freq_ind) * 1.05]; % ±5% bandwidth
                [b, a] = butter(2, [min(fb)/Fn max(fb)/Fn]); % Butterworth filter design
                filtered_signal = filtfilt(b, a, Vm_trace); % filtfilt
                analytic_signal = hilbert(filtered_signal); % Compute the analytic signal
                wavelet_matrix(freq_ind, :, trial_ind) = angle(analytic_signal); 
            end
        end
        
        for freq_ind = 1:length(frequencies)
            for time_ind = 1:size(wavelet_matrix, 2)
                phase_data = exp(1i * squeeze(wavelet_matrix(freq_ind, time_ind, :))); % Extract phase data
                itc_value = abs(mean(phase_data, 'omitnan')); % Compute ITC as the mean vector length
                itc_values(freq_ind, time_ind) = itc_values(freq_ind, time_ind) + itc_value; % Sum across neurons
            end
        end
    end
    itc_values = itc_values / length(neurons); % Average across neurons
end


function [itc_mean, itc_sem] = compute_itc_sem(neurons, vm_all, frequencies, num_trials, time_window)
    FS = 1000;
    Fn = FS / 2;
    itc_all = zeros(length(neurons), length(frequencies), FS * time_window + 1);
    
    % Loop through neurons
    for cell_ind = 1:length(neurons)
        curr_cell = neurons(cell_ind); % Current neuron index
        Vm_curr_cell = squeeze(vm_all(:, :, curr_cell)); % Extract data for the current neuron
        Vm_curr_cell(isnan(Vm_curr_cell)) = 0; % Replace NaNs with zero

        wavelet_matrix = zeros(length(frequencies), FS * time_window + 1, num_trials(cell_ind));

        % Loop through trials
        for trial_ind = 1:num_trials(cell_ind)
            Vm_trace = Vm_curr_cell(trial_ind, :); % Extract voltage trace for this trial
            Vm_trace(isnan(Vm_trace)) = 0; % Replace NaNs with zero
            
            % Loop through frequencies
            for freq_ind = 1:length(frequencies)
                fb = [frequencies(freq_ind) * 0.95, frequencies(freq_ind) * 1.05]; % ±5% bandwidth
                [b, a] = butter(2, [min(fb)/Fn max(fb)/Fn]); % Butterworth filter design
                filtered_signal = filtfilt(b, a, Vm_trace); % Apply the filter
                analytic_signal = hilbert(filtered_signal); % Compute the analytic signal
                wavelet_matrix(freq_ind, :, trial_ind) = angle(analytic_signal); % Store phase
            end
        end

        % Compute ITC for the current neuron
        for freq_ind = 1:length(frequencies)
            for time_ind = 1:size(wavelet_matrix, 2)
                phase_data = exp(1i * squeeze(wavelet_matrix(freq_ind, time_ind, :))); % Extract phase data
                itc_value = abs(mean(phase_data, 'omitnan')); % Compute ITC as the mean vector length
                itc_all(cell_ind, freq_ind, time_ind) = itc_value; % Store ITC for this neuron
            end
        end
    end

    % Compute mean and SEM across neurons
    itc_mean = squeeze(mean(itc_all, 1, 'omitnan')); % Mean ITC across neurons
    itc_sem = squeeze(std(itc_all, 0, 1, 'omitnan') ./ sqrt(size(itc_all, 1))); % SEM across neurons
end


