%This script file uses the Fieldtrip toolbox to compute night-to-night
%REM (Theeta) spectral analysis results for the RESTORE dataset.

clear
clc
close all

%% User settings
sqi_cutoff = 3; %Rejects epochs if the log-standard deviation of the EEG in that epoch is more than sqi_cutoff standard deviations from the mean log-standard deviation per-epoch for that sleep state

%Paths
fieldtrip_path = './fieldtrip'; %Path to fieldtrip toolbox
record_path = './restore'; %Path to RESTORE files
record_path_suffix = {'Bsl', 'Int_1', 'Int_2', 'Int_3', 'Int_4', 'Int_5', 'Int_6', 'Int_7', 'Flw'}; %Suffixes used to specify protocol night

%Loads the FieldTrip Toolbox
addpath(fieldtrip_path);

%% Pre-allocates arrays and settings
bpw_ngt_avg = nan(9, 19); %Overnight bandpower, average
bpw_ngt_std = nan(9, 19); %Overnight bandpower, standard deviation
npw_ngt_avg = nan(9, 19); %Percentage overnight bandpower, average
npw_ngt_std = nan(9, 19); %Percentage overnight bandpower, standard deviation

%Filter settings
fs = 256;
[b, a] = butter(5, [0.3/(fs/2) 35.0/(fs/2)]);

%% Iterates through each night of SRT
for i = 1:9

    %Gets the list of participants with data for the selected night
    records_night = dir([record_path '/' record_path_suffix{i}]);
    records_night = records_night(~startsWith({records_night.name}, '.'));
    
    %Iterates through each participant's data
    for j = 1:length(records_night)
        disp(['Now processing record ' records_night(j).name])
        part_num = str2double(records_night(j).name(1:3));

        %Loads the file and its annotations
        edf_meta = edfinfo([records_night(j).folder '/' records_night(j).name]);
        edf_file = edfread([records_night(j).folder '/' records_night(j).name]);

        %Loads annotation file
        anno = readtable([record_path '/hypnograms/' records_night(j).name(1:20) '_' record_path_suffix{i}(1:3) '.csv']);
        trim = readmatrix([record_path '/trimming_rem/'  records_night(j).name(1:20) '_' record_path_suffix{i}(1:3) '_scoring_REM.txt']);
        if trim(1, 1) == 6 %Skips artefact epochs at the beginning of the hypnogram
            trim_start = find(trim(:,1) ~= 6, 1, 'first');
            trim = trim(trim_start:end, :);
        end

        %% Processes Annotations

        %Reformats manual annotations
        stages = -ones(height(anno), 2);
        stages(:,1) = seconds(anno.SignalID - anno.SignalID(1));
        if length(stages) == length(trim)
            stages(:,2) = trim(:, 1);
        else
            disp('ERROR: Trimming array inconsistent with hypnogram')
            return
        end
    
        %% Prepares the EEG signal for analysis
        if length(edf_file.Fz{1}) ~= fs
            disp('ERROR: Incorrect Sample Rate Specified')
            return
        end
        eeg = cell2mat(edf_file.Fz);

        % Deals with offsets
        eeg_start = datetime(strcat(edf_meta.StartDate, ".", edf_meta.StartTime), 'InputFormat', 'dd.MM.yy.HH.mm.ss');
        off_str = eeg_start - anno.SignalID(1); %eeg starts off_str seconds later than annotation
        
        %Adjusts for cases where the EEG begins after annotation
        if seconds(off_str) > 0
            bad_rows = ceil(seconds(off_str)/30);
            stages(1:bad_rows, :) = [];
            anno(1:bad_rows, :) = [];
            stages(: ,1) = stages(:, 1)-30*bad_rows;
            off_str = eeg_start - anno.SignalID(1);
        end
        
        %Crops the EEG signal to be the same length as the annotations
        off_end = (eeg_start + seconds(length(eeg)/fs)) - (anno.SignalID(end) + seconds(30)); %eeg ends off_end seconds later than annotation
        if seconds(off_str) <= 0
            eeg(1:fs*seconds(-off_str)) = [];
            eeg(end-fs*seconds(off_end)+1:end) = [];
        else
            disp('ERROR: EEG recording begins after annotation')
            return
        end

        %Checks that signal lengths are now consistent
        if length(eeg)/fs ~= length(stages)*30
            disp('ERROR: EEG and annotation lengths are not consistent after clipping')
            return
        end

        %% Performs Standard Deviation Based Signal Quality Check
        stages = stages(trim(:,2) == 0, :); %Applies manual trimming to ensure comparable amounts of NREM sleep
        epoch_lbls = stages(:,2);

        %Gets standard deviation for each epoch of NREM sleep
        epoch_std = zeros(length(epoch_lbls), 1);
        for k = 1:length(epoch_lbls)
            if epoch_lbls(k) == 5
                st = (k - 1)*fs*30 + 1;
                et = (k    )*fs*30;
                if et <= length(eeg)
                    eeg_seg = eeg(st:et);
                    epoch_std(k) = log(std(eeg_seg));
                else
                    disp('ERROR: EEG file length error')
                    return
                end
            end
        end

        %Computes the distribution of standard deviations for each group of NREM sleep
        epoch_sqi = (sqi_cutoff+1)*ones(length(epoch_lbls), 1);
        for k = 5
            stage_loc = find(epoch_lbls == k);
            stage_std = zscore(epoch_std(epoch_lbls == k));
            epoch_sqi(stage_loc) = stage_std;
        end

        %Uses the SQI to label noisy epochs
        stages(epoch_sqi >= sqi_cutoff, 2) = -1;
        disp(['Rejected ' num2str(sum(epoch_sqi >= sqi_cutoff)) ' epochs based on log standard deviation'])

        %% Gets a list of contiguous REM segment start and end times
        stages = stages(stages(:,2) > 0, :);
        stage_diffs = diff(stages(:, 1));
        epoch_gaps = [true; stage_diffs > 30];
        rem_segment_ids = cumsum(epoch_gaps);
        rem_segments = splitapply(@(t) [min(t), max(t) + 30], stages(:, 1), rem_segment_ids);

        %% Uses Fieldtrip to Perform Spectral Analysis
        eeg = filtfilt(b, a, eeg)';

        %Creates a Fieldtrip object containing the EEG signal
        eeg_data = [];
        eeg_data.label{1} = 'Fz';
        eeg_data.time{1} = (0:1/fs:(length(eeg)-1)/fs);
        eeg_data.trial{1} = eeg;
        eeg_data.fsample = fs;
        eeg_data.sampleinfo = [1 length(eeg)];

        %Breaks the EEG signal into contiguous NREM segments
        eeg_segs = struct();
        eeg_segs.trl = [rem_segments(:, 1)*fs+1, rem_segments(:, 2)*fs, zeros(size(rem_segments, 1), 1)];
        eeg_data = ft_redefinetrial(eeg_segs, eeg_data);

        %Uses Fieldtrip to split the EEG into 4 second windows, with a 50% overlap
        cfg = [];
        cfg.length = 4;
        cfg.overlap = 0.5;
        eeg_data = ft_redefinetrial(cfg, eeg_data);

        %Fieldtrip configuration for spectral analysis
        cfg_spectral = [];
        cfg_spectral.method     = 'mtmfft';
        cfg_spectral.output     = 'pow';
        cfg_spectral.taper      = 'hanning';  
        cfg_spectral.pad        = 'nextpow2'; 
        cfg_spectral.foilim     = [0.5 32];    
        cfg_spectral.keeptrials = 'yes';

        %Spectral analysis using fieldtrip
        freq = ft_freqanalysis(cfg_spectral, eeg_data);
        
        %Computes bandpower and normalised bandpower across the night
        bpw_epoch = zeros(length(freq.cumsumcnt), 1);
        npw_epoch = zeros(length(freq.cumsumcnt), 1);
        for k = 1:length(freq.cumsumcnt)
            bpw_epoch(k) = simps(freq.freq( 17: 31), freq.powspctrm(k, :,  17: 31)); %theta
            tot_power = simps(freq.freq(1:127), freq.powspctrm(k, :, 1:127, :));
            npw_epoch(k) = bpw_epoch(k)./tot_power*100;
        end

        %% Aggregates Overall PSD and BandPower
        bpw_ngt_avg(i, part_num) = mean(bpw_epoch, 1);
        bpw_ngt_std(i, part_num) = std(bpw_epoch, [], 1);
        npw_ngt_avg(i, part_num) = mean(npw_epoch, 1);
        npw_ngt_std(i, part_num) = std(npw_epoch, [], 1);
    end
end

%% Prepares Results for Plotting and Writing to Disk
bpw_ovr_avg = nanmean(bpw_ngt_avg, 2);
bpw_ovr_std  = nanstd(bpw_ngt_avg, [], 2);
npw_ovr_avg = nanmean(npw_ngt_avg, 2);
npw_ovr_std  = nanstd(npw_ngt_avg, [], 2);

%% Plots distribution of band power for all participants
band_labels  = {'Theta'};
night_labels = {'B1', 'I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'F1'};

figure(1)
set(gcf, 'units', 'centimeters', 'Position', [15 15 15.5 11.625])
tiledlayout(1, 1,'TileSpacing','Compact','Padding','Compact');

for i = 1
    nexttile
    set(gca, 'FontName', 'Times New Roman')
    hold on
    title(band_labels(i))
    plot(ones(19, 1)*(1:9), squeeze(bpw_ngt_avg(:, :, i))', 'kx')
    errorbar(1:9, bpw_ovr_avg(:, i), bpw_ovr_std(:, i), 'b-');
    xticks(1:9)
    xlim([0 10])
    xticklabels(night_labels)
    xlabel('Study Night')
    ylabel('Band Power, μV^{2}')
end

%% Plots distribution of percentage power for all participants
figure(2)
set(gcf, 'units', 'centimeters', 'Position', [15 15 15.5 11.625])
tiledlayout(1, 1,'TileSpacing','Compact','Padding','Compact');

for i = 1
    nexttile
    set(gca, 'FontName', 'Times New Roman')
    hold on
    title(band_labels(i))
    plot(ones(19, 1)*(1:9), squeeze(npw_ngt_avg(:, :, i))', 'kx')
    errorbar(1:9, npw_ovr_avg(:, i), npw_ovr_std(:, i), 'b-');
    xticks(1:9)
    xlim([0 10])
    xticklabels(night_labels)
    xlabel('Study Night')
    ylabel('Percentage Power (%)')
end

%% Saves figures
w = 15.5;
h = 17.4375;

figure(1)
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'centimeters');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- w)/2;
bottom = (papersize(2)- h)/2;
myfiguresize = [left, bottom, w, h];
set(gcf,'PaperPosition', myfiguresize);
pause(1)
print('./figures/restore/rest_bpw_ft_rem', '-dpng', '-r300');

figure(2)
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'centimeters');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- w)/2;
bottom = (papersize(2)- h)/2;
myfiguresize = [left, bottom, w, h];
set(gcf,'PaperPosition', myfiguresize);
pause(1)
print('./figures/restore/rest_npw_ft_rem', '-dpng', '-r300');

%% Writes Spectral Results to Disk
save('./figures/restore/rest_spectral_ft_rem', 'npw_ngt_avg', 'bpw_ngt_avg', 'npw_ngt_std', 'bpw_ngt_std', 'npw_ovr_avg', 'bpw_ovr_avg', 'npw_ovr_std', 'bpw_ovr_std')