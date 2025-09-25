%% slopes for responsive neurons (all) 

FS = 1000;
numTrials = 40;

% Define analysis windows
baseline_idx  = 501:1000;
transient_idx = 1001:1250;
post_idx      = 2001:3000;

% Select only transiently modulated neurons
numNeurons = length(good_neurons);

z_baseline  = NaN(numTrials, numNeurons);
z_transient = NaN(numTrials, numNeurons);
z_post      = NaN(numTrials, numNeurons);

% Z-score within each trial
for n = 1:numNeurons
    neuronID = n;
    for t = 1:numTrials
        base_data = vm_all(t, 1:1000, neuronID);
        base_mean = mean(base_data, 'omitnan');
        base_std  = std(base_data, 'omitnan');
        if base_std == 0
            base_std = 1;
        end

        z_baseline(t, n)  = (mean(vm_all(t, baseline_idx,  neuronID), 'omitnan') - base_mean) / base_std;
        z_transient(t, n) = (mean(vm_all(t, transient_idx, neuronID), 'omitnan') - base_mean) / base_std;
        z_post(t, n)      = (mean(vm_all(t, post_idx,      neuronID), 'omitnan') - base_mean) / base_std;
    end
end

% Compute slopes across trials
x = (1:numTrials)';
slopes_baseline  = NaN(1, numNeurons);
slopes_transient = NaN(1, numNeurons);
slopes_post      = NaN(1, numNeurons);

for n = 1:numNeurons
    % Baseline
    y_b = z_baseline(:, n);
    if sum(~isnan(y_b)) > 1
        p_b = polyfit(x(~isnan(y_b)), y_b(~isnan(y_b)), 1);
        slopes_baseline(n) = p_b(1);
    end

    % Transient
    y_t = z_transient(:, n);
    if sum(~isnan(y_t)) > 1
        p_t = polyfit(x(~isnan(y_t)), y_t(~isnan(y_t)), 1);
        slopes_transient(n) = p_t(1);
    end

    % Post
    y_p = z_post(:, n);
    if sum(~isnan(y_p)) > 1
        p_p = polyfit(x(~isnan(y_p)), y_p(~isnan(y_p)), 1);
        slopes_post(n) = p_p(1);
    end
end

% Shared histogram bin edges
all_slopes = [slopes_baseline, slopes_transient, slopes_post];
bin_edges = linspace(min(all_slopes), max(all_slopes), 30);

% Plot histograms
figure;

subplot(3,1,1);
histogram(slopes_baseline, 'BinEdges', bin_edges, 'FaceColor', 'k');
xlabel('Slope of Z-scored Response');
ylabel('Number of Neurons');
title('Baseline Slopes');
grid on;
%ylim([0 60]);

subplot(3,1,2);
histogram(slopes_transient, 'BinEdges', bin_edges, 'FaceColor', 'b');
xlabel('Slope of Z-scored Response');
ylabel('Number of Neurons');
title('Transient Slopes');
grid on;
%ylim([0 25]);

subplot(3,1,3);
histogram(slopes_post, 'BinEdges', bin_edges, 'FaceColor', 'r');
xlabel('Slope of Z-scored Response');
ylabel('Number of Neurons');
title('Post Slopes');
grid on;
%ylim([0 60]);

%%

[p_signedrank, ~, stats_signedrank] = signrank(slopes.transient);
fprintf('Wilcoxon signed-rank test: W = %d, p = %.4f\n', stats_signedrank.signedrank, p_signedrank);

[p_signtest, ~] = signtest(slopes.transient, 0);
fprintf('Sign test: p = %.4f\n', p_signtest);

skew_val = skewness(slopes.transient);
fprintf('Skewness of transient slopes: %.4f\n', skew_val);

[h_norm, p_norm] = kstest((slopes.transient - mean(slopes.transient)) / std(slopes.transient));
fprintf('Kolmogorov–Smirnov test for normality: p = %.4f\n', p_norm);

%Interpretation:
%- The median transient slope is not significantly different from zero (p = 0.2107).
%- However, the skewness (-0.3895) suggests a mild trend toward more decreasing responses.

%% specific neuron response across trials (neuron 50)
FS = 1000;
numTrials = 20;
neuronID = 39; % 39 and 112

vm_trials = squeeze(vm_all(:, :, 10));

valid_trials = all(~isnan(vm_trials(:, 100)), 2);

baseline_idx        = 501:1000;
transient_idx       = 1001:1250;
delayed_idx         = 1251:2000;
post_transient_idx  = 2001:2250;
post_delayed_idx    = 2251:3000;

resp_baseline       = NaN(numTrials, 1);
resp_transient      = NaN(numTrials, 1);
resp_delayed        = NaN(numTrials, 1);
resp_post_transient = NaN(numTrials, 1);
resp_post_delayed   = NaN(numTrials, 1);

% Z-score and extract trial-by-trial responses
for t = 1:numTrials
    baseline_data = vm_all(t, 1:1000, neuronID);
    base_mean = mean(baseline_data, 'omitnan');
    base_std  = std(baseline_data,  'omitnan');
    if base_std == 0
        base_std = 1;
    end

    resp_baseline(t)       = (mean(vm_all(t, baseline_idx,        neuronID), 'omitnan') - base_mean) / base_std;
    resp_transient(t)      = (mean(vm_all(t, transient_idx,       neuronID), 'omitnan') - base_mean) / base_std;
    resp_delayed(t)        = (mean(vm_all(t, delayed_idx,         neuronID), 'omitnan') - base_mean) / base_std;
    resp_post_transient(t) = (mean(vm_all(t, post_transient_idx,  neuronID), 'omitnan') - base_mean) / base_std;
    resp_post_delayed(t)   = (mean(vm_all(t, post_delayed_idx,    neuronID), 'omitnan') - base_mean) / base_std;
end

% Compute linear fits
x = (1:numTrials)';
p_trans      = polyfit(x, resp_transient, 1);
p_delayed    = polyfit(x, resp_delayed, 1);
p_post_trans = polyfit(x, resp_post_transient, 1);
p_post_delay = polyfit(x, resp_post_delayed, 1);

% Predicted fit lines
fit_trans      = polyval(p_trans, x);
fit_delayed    = polyval(p_delayed, x);
fit_post_trans = polyval(p_post_trans, x);
fit_post_delay = polyval(p_post_delay, x);

figure; hold on;

h1 = plot(x, resp_transient,      '-b', 'LineWidth', 1.5);
h2 = plot(x, fit_trans,           '--b', 'LineWidth', 1.5);

h3 = plot(x, resp_post_delayed,        '-g', 'LineWidth', 1.5);
h4 = plot(x, fit_post_delay,         '--g', 'LineWidth', 1.5);

xlabel('Trial Number');
ylabel('Z-scored Vm Response');
title(sprintf('Trial-by-Trial Responses (Neuron %d)', neuronID));

legend([h1, h2, h3, h4], {
    sprintf('Transient (slope = %.3f)',      p_trans(1)), ...
    'Fit: Transient', ...
    sprintf('Delayed (slope = %.3f)',        p_delayed(1)), ...
    'Fit: Delayed', ...
}, 'Location', 'best');

