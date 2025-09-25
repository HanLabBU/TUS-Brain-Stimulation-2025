%% Plot all neurons raster
for i = 1:161
    
    curr_mod_roi = i;

    neuronPlot = figure(1);
    
    plot_spike_raster(neuronPlot,full_block_trace,full_block_raster,full_block_Vm,full_block_psth,preIdx,curr_mod_roi,T_params, mod_Sp_track_250ms, mod_Sp_track_1000ms, mod_Vm_track_250ms, mod_Vm_track_1000ms)
    
    fprintf('\nPress any key to continue')
    
    pause
end

%% function
function [] = plot_spike_raster(handle_fig, full_block_trace,full_block_raster,full_block_Vm,full_block_psth,preIdx,curr_roi,T_params, mod_Sp_track_250ms, mod_Sp_track_1000ms, mod_Vm_track_250ms, mod_Vm_track_1000ms)
    figure(handle_fig) 
    FS = 1000;
    sgfilt_wind = 21;
    
    Mod_Status_Curr_ROI = [...
        ', Sp Mod, 250ms: ',num2str(mod_Sp_track_250ms(curr_roi)),...
        ', Sp Mod, 1000ms: ',num2str(mod_Sp_track_1000ms(curr_roi)),...
        ', Vm Mod, 250ms: ',num2str(mod_Vm_track_250ms(curr_roi)),...
        ', Vm Mod, 1000ms: ',num2str(mod_Vm_track_1000ms(curr_roi))];

    aa = tiledlayout(4,1);
    title(aa,Mod_Status_Curr_ROI)

    nexttile(2,[2,1]);
    hold on; 
    
    trial_curr = 0;
    psth = squeeze(full_block_psth(:,:,curr_roi)); 
    x_ax = ((1:size(full_block_trace,2))-1)./FS-(preIdx+1)/FS;

    for tt = 1:size(full_block_trace,1) % loop through trials 
        if sum(isnan(full_block_trace(tt,:,curr_roi))) == size(full_block_trace,2)
            continue % trial isn't good
        end
        
        trial_curr = trial_curr+1;
        spikes = squeeze(full_block_raster(tt,:,curr_roi));    
        spikes(isnan(spikes)) = 0;
        spike_inds = find(spikes);    
        trace_detrend_norm = squeeze(full_block_trace(tt,:,curr_roi));

        trace_Vm = squeeze(full_block_Vm(tt,:,curr_roi));
    
        plot(x_ax,trace_detrend_norm + 1*(trial_curr-1),'k')   

        trace_amp = max(trace_detrend_norm,[],'omitnan');
        plot(x_ax(spike_inds),trace_amp*ones(size(spike_inds)) + 1*(trial_curr-1),'.','Color',[255 36 0]./255)
    end    

    psth(:,1:30) = NaN;
    psth(:,end-30:end) = NaN;    

    % 95 percent CI 
    psth_high_CI_bound = mean(psth,1,'omitnan')+1.96*std(psth,[],1,'omitnan')./sqrt(size(psth,1));
    psth_low_CI_bound = mean(psth,1,'omitnan')-1.96*std(psth,[],1,'omitnan')./sqrt(size(psth,1));
    psth_CI_fill = [psth_high_CI_bound'; flipud(psth_low_CI_bound')]';
    t2 = [x_ax,fliplr(x_ax)];
    t2(isnan(psth_CI_fill)) = [];
    psth_CI_fill(isnan(psth_CI_fill)) = [];   

    xline([0 1],'color',[0.3 0.6 1])
    
    title(sprintf('ROI %i',curr_roi))
    xlabel('Time (s)')
    yticks((0:1:(1*(tt-1))))
    yticklabels(1:tt)
    ylabel('Trial')
    axis tight
    xlim([-1 2])

    nexttile(1,[1,1]);
    hold on; 
    fill(t2,psth_CI_fill,[255 36 0]./255,'EdgeColor','none','FaceAlpha',0.15)
    plot(x_ax,mean(psth,1,'omitnan'),'Color',[255 36 0]./255)
    xline([0 1],'color',[0.3 0.6 1])
    
    title(sprintf('ROI %i',curr_roi))
    xlabel('Time (s)')
    ylabel('Spike Rate (Hz)')
    axis tight
    xlim([-1 2])

    % 95 percent CI 
    Vm = squeeze(full_block_Vm(:,:,curr_roi));
    Vm_high_CI_bound = mean(Vm,1,'omitnan')+1.96*std(Vm,[],1,'omitnan')./sqrt(size(Vm,1));
    Vm_low_CI_bound = mean(Vm,1,'omitnan')-1.96*std(Vm,[],1,'omitnan')./sqrt(size(Vm,1));
    Vm_CI_fill = [Vm_high_CI_bound'; flipud(Vm_low_CI_bound')]';
    t2 = [x_ax,fliplr(x_ax)];
    t2(isnan(Vm_CI_fill)) = [];
    Vm_CI_fill(isnan(Vm_CI_fill)) = [];   

    nexttile(4,[1,1]);
    hold on; 
    mean_Vm = mean(Vm,1,'omitnan');
    fill(t2,Vm_CI_fill,[255 36 0]./255,'EdgeColor','none','FaceAlpha',0.15)
    plot(x_ax,mean_Vm,'Color','r')
    xline([0 1],'color',[0.3 0.6 1])
    
    title(sprintf('ROI %i',curr_roi))
    xlabel('Time (s)')
    ylabel('Avg. Vm')
    axis tight
    xlim([-1 2])
end