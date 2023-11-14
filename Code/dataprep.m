%% Code Outline

% This code does the following:
%       - loads in both fMRI and behavioural data
%       - processing both behavioural and fMRI data
%       - removes subjects that don't meet selection criteria
%       - normalises amount of trials for each subject
%       - checks robustness of random sampling


%% Load in task data

% Baseline + Rotation
rootdir = 'D:\PhD\visuomotor_Gallivan\BIDS'; %set path of root directory
filelist = dir(fullfile(rootdir, '**\*ses-01_task-rotation_events.tsv*')); %get list of .tsv files for all subfolders

% Create structure to hold the data
storage(length(filelist)) = struct('name',1,'duration',1,'onset',1,'response_error',1,'response_time',1,'trial_number',1,'trial_type',1);
% Load in each subject and store their data
for ii = 1:length(filelist)
    subjectname = extractBefore(filelist(ii).name,'_task-rotation_events.tsv'); 
    storage(ii).name = subjectname;
    fullFileName = fullfile(filelist(ii).folder, filelist(ii).name); %create absolute path to tsv file
    tdfread(fullFileName); %read in tsv file
    storage(ii).duration = duration;
    storage(ii).onset = onset;
    storage(ii).response_error = response_error;
    storage(ii).response_time = response_time;
    storage(ii).trial_number = trial_number;
    storage(ii).trial_type = trial_type;
end

% Washout
filelist2 = dir(fullfile(rootdir, '**\*ses-01_task-washout_events.tsv*'));
for ii = 1:length(filelist2)
    fullFileName = fullfile(filelist2(ii).folder, filelist2(ii).name);
    tdfread(fullFileName);
    storage(ii).duration2 = duration;
    storage(ii).onset2 = onset;
    storage(ii).response_error2 = response_error;
    storage(ii).response_time2 = response_time;
    storage(ii).trial_number2 = trial_number;
    storage(ii).trial_type2 = trial_type;
end


%% Load in time-series data for Day 1

% Baseline + Rotation
% Specify the folder where the files lives.
myFolder = 'D:\PhD\visuomotor_Gallivan\denoise';

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, 'voltron_400\baseline_rotation\*ses-01*.mat');
theFiles = dir(filePattern); % Storing all filenames

% Loop to load in the data, store name and values of zscore
for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name; 
    fullFileName = fullfile(theFiles(ii).folder, baseFileName); %making absolute path to file
    ts = load(fullFileName); %loads in data
    storage(ii).ts_baserot = ts.ts; %stores data under .data
end

% Washout
% Get a list of all files in the folder with the desired file name pattern.
filePattern2 = fullfile(myFolder, 'voltron_400\washout\*ses-01*.mat');
theFiles2 = dir(filePattern2); % Storing all filenames

% Loop to load in the data, store name and values of zscore
for ii = 1:length(theFiles2)
    baseFileName = theFiles2(ii).name; 
    fullFileName = fullfile(theFiles2(ii).folder, baseFileName); %making absolute path to file
    ts = load(fullFileName); %loads in data
    storage(ii).ts_wash = ts.ts; %stores data under .data
end


%% Data processing

% response error can have n/a
% response time same behaviour as above
% 17 subjects have n/a in responses

% Some subjects have trials with 0 response time
% RT removed: RT>2s or RT<100ms
% RE removed: RE>100

% Baseline + Rotation: Change response_time to double & change n/a to NaN
for ii = 1:length(storage)
    RT = storage(ii).response_time;
    double_test = zeros(length(RT),1);
    if isa(RT,'double')
        storage(ii).RT = RT;
    else
        for jj = 1:length(RT)
            num = str2double(RT(jj,:));
            double_test(jj) = num;
        end
        storage(ii).RT = double_test;
    end
end

% Washout: Change response_time to double & change n/a to NaN
for ii = 1:length(storage)
    RT = storage(ii).response_time2;
    double_test = zeros(length(RT),1);
    if isa(RT,'double')
        storage(ii).RT2 = RT;
    else
        for jj = 1:length(RT)
            num = str2double(RT(jj,:));
            double_test(jj) = num;
        end
        storage(ii).RT2 = double_test;
    end
end

% Baseline + Rotation: Change response_error to double & change n/a to NaN
for ii = 1:length(storage)
    RE = storage(ii).response_error;
    double_test = zeros(length(RE),1);
    if isa(RE,'double')
        storage(ii).RE = RE;
    else
        for jj = 1:length(RE)
            num = str2double(RE(jj,:));
            double_test(jj) = num;
        end
        storage(ii).RE = double_test;
    end
end

% Washout: Change response_error to double & change n/a to NaN
for ii = 1:length(storage)
    RE = storage(ii).response_error2;
    double_test = zeros(length(RE),1);
    if isa(RE,'double')
        storage(ii).RE2 = RE;
    else
        for jj = 1:length(RE)
            num = str2double(RE(jj,:));
            double_test(jj) = num;
        end
        storage(ii).RE2 = double_test;
    end
end

% Add set number for each trial
set = 1:55;
set2 = repelem(set,8)';
for ii = 1:length(storage)
    trial = storage(ii).trial_number;
    trial(:,2) = set2;
    storage(ii).trial_number = trial;
end

set = 56:70;
set2 = repelem(set,8)';
for ii = 1:length(storage)
    trial = storage(ii).trial_number2;
    trial(:,2) = set2;
    storage(ii).trial_number2 = trial;
end


%% Changing outliers to NaN

% Baseline + Rotation
% Create trial_type factor variable (same across subjects)
trial_type = storage(1).trial_type;
trial_type_f = zeros(440,1);
for ii = 1:length(trial_type)
    if isequal(trial_type(ii,:),'baseline')
        trial_type_f(ii) = 0;
    elseif isequal(trial_type(ii,:),'rotation')
        trial_type_f(ii) = 1;
    end
end

% Change outliers to NaN (based off RT)
for ii = 1:length(storage)
    task_subject = zeros(440,5);
    task_subject(:,1) = storage(ii).trial_number(:,1);
    task_subject(:,2) = trial_type_f;
    RT = storage(ii).RT;
    RE = abs(storage(ii).RE);
    new_RT = RT;
    new_RT(RT>2 | RT<0.3) = NaN;
    new_RE = RE;
    new_RE(RE>100) = NaN;
    task_subject(:,3) = new_RT;
    task_subject(:,4) = new_RE;
    task_subject(:,5) = storage(ii).trial_number(:,2);
    storage(ii).dm = task_subject;
end

% Washout
% Create trial_type factor variable (same across subjects)
trial_type2 = storage(1).trial_type2;
trial_type_f2 = zeros(120,1);
for ii = 1:length(trial_type2)
    if isequal(trial_type2(ii,:),'washout')
        trial_type_f2(ii) = 2;
    end
end

% Change outliers to NaN (based off RT)
for ii = 1:length(storage)
    task_subject = zeros(120,5);
    task_subject(:,1) = storage(ii).trial_number2(:,1);
    task_subject(:,2) = trial_type_f2;
    RT = storage(ii).RT2;
    RE = abs(storage(ii).RE2);
    new_RT = RT;
    new_RT(RT>2 | RT<0.3) = NaN;
    new_RE = RE;
    new_RE(RE>100) = NaN;
    task_subject(:,3) = new_RT;
    task_subject(:,4) = new_RE;
    task_subject(:,5) = storage(ii).trial_number2(:,2);
    storage(ii).dm_wash = task_subject;
end

% Combine dm
for ii = 1:length(storage)
    dm = storage(ii).dm;
    dm_wash = storage(ii).dm_wash;
    storage(ii).dm_all = [dm; dm_wash];
end


%% Extract TRs for cortical and cerebellar ROIs (ignoring missing and outlier TRs)

% Baseline + Rotation
% Cortex, cerebellar, basal ganglia time-series (with FILL TRs)
% Basal ganglia = dorsal striatum (caudate nucleus, putamen), ventral
% striatum (nucleus accumbens, olfactory tubercle), globus pallidus,
% ventral pallidum, substantia nigra, and subthalamic nuclei
for ii = 1:length(storage)
    ts = storage(ii).ts_baserot;
    ts_cortex = ts(:,1:400);
    ts_cerebellum = ts(:,455:482);
    ts_cnR = ts(:,413:420);
    ts_cnL = ts(:,440:447);
    ts_nacR = ts(:,424:425);
    ts_nacL = ts(:,451:452);
    ts_gpR = ts(:,426:427);
    ts_gpL = ts(:,453:454);
    ts_stn = ts(:,499:500);
    ts_bg = [ts_cnR ts_nacR ts_gpR ts_cnL ts_nacL ts_gpL ts_stn];
    storage(ii).ts_baserot_roi = [ts_cortex ts_cerebellum ts_bg];
end

% Washout
% Cortex, cerebellar, basal ganglia time-series (with FILL TRs)
for ii = 1:length(storage)
    ts = storage(ii).ts_wash;
    ts_cortex = ts(:,1:400);
    ts_cerebellum = ts(:,455:482);
    ts_cnR = ts(:,413:420);
    ts_cnL = ts(:,440:447);
    ts_nacR = ts(:,424:425);
    ts_nacL = ts(:,451:452);
    ts_gpR = ts(:,426:427);
    ts_gpL = ts(:,453:454);
    ts_stn = ts(:,499:500);
    ts_bg = [ts_cnR ts_nacR ts_gpR ts_cnL ts_nacL ts_gpL ts_stn];
    storage(ii).ts_wash_roi = [ts_cortex ts_cerebellum ts_bg];
end


%% Explore dataset

% Clear up storage
storage = rmfield(storage, 'duration');
storage = rmfield(storage, 'duration2');
storage = rmfield(storage, 'onset');
storage = rmfield(storage, 'onset2');
storage = rmfield(storage, 'response_error2');
storage = rmfield(storage, 'response_error');
storage = rmfield(storage, 'response_time2');
storage = rmfield(storage, 'response_time');
storage = rmfield(storage, 'trial_number');
storage = rmfield(storage, 'trial_number2');
storage = rmfield(storage, 'trial_type');
storage = rmfield(storage, 'trial_type2');
storage = rmfield(storage, 'RT');
storage = rmfield(storage, 'RT2');
storage = rmfield(storage, 'RE');
storage = rmfield(storage, 'RE2');
storage = rmfield(storage, 'ts_baserot');
storage = rmfield(storage, 'ts_wash');

% Plot RT across trials
figure;
for ii = 1:length(storage)
    dm = storage(ii).dm;
    RT = dm(:,3);
    scatter(1:length(RT),RT,30,'filled');
    xlabel('Trial')
    ylabel('RT')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    xline(120)
    xline(440)
    pause
end

% Plot RT across sets
all_set = storage(1).dm(:,5);
set_name = unique(all_set);
figure;
for ii = 1:length(storage)
    dm = storage(ii).dm;
    RT = dm(:,3);
    for jj = 1:length(set_name)
        RT_set(jj) = mean(RT(set_name==jj));
    end
    storage(ii).RT_set = RT_set;
    scatter(1:length(set_name),RT_set,30,'filled');
    xlabel('Set')
    ylabel('RT')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    xline(15)
    xline(55)
    pause
end

% Plot RE across trials
figure;
for ii = 1:length(storage)
    dm = storage(ii).dm;
    RE = dm(:,4);
    scatter(1:length(RE),RE,30,'filled');
    xlabel('Trial')
    ylabel('RE')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    xline(120)
    xline(440)
    pause
end

% Plot RE across sets
figure;
for ii = 1:length(storage)
    dm = storage(ii).dm;
    RE = dm(:,4);
    for jj = 1:length(set_name)
        RE_set(jj) = mean(RE(set_name==jj));
    end
    storage(ii).RE_set = RE_set;
    scatter(1:length(set_name),RE_set,30,'filled');
    xlabel('Set')
    ylabel('RE')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    xline(15)
    xline(55)
    pause
end

% Plot RT vs. RE for baseline
all_dm = vertcat(storage.dm);
min_RT = min(all_dm(:,3));
max_RT = max(all_dm(:,3));
min_RE = min(all_dm(:,4));
max_RE = max(all_dm(:,4));
base = 120;
figure;
for ii = 1:length(storage)
    dm = storage(ii).dm;
    RT = dm(1:base,3);
    RE = dm(1:base,4);
    scatter(RT,RE,30,1:length(RT),'filled');
    colormap(flipud(hot));
    xlabel('RT')
    ylabel('RE')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    xlim([min_RT max_RT]);
    ylim([min_RE max_RE]);
    pause
end

% Plot RT vs. RE for rotation
rot = 440;
figure;
for ii = 1:length(storage)
    dm = storage(ii).dm;
    RT = dm(base+1:rot,3);
    RE = dm(base+1:rot,4);
    scatter(RT,RE,30,1:length(RT),'filled');
    colormap(flipud(hot));
    xlabel('RT')
    ylabel('RE')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    xlim([min_RT max_RT]);
    ylim([min_RE max_RE]);
    pause
end

% Plot RT vs. RE for washout
wash = 560;
figure;
for ii = 1:length(storage)
    dm = storage(ii).dm;
    RT = dm(rot+1:wash,3);
    RE = dm(rot+1:wash,4);
    scatter(RT,RE,30,1:length(RT),'filled');
    colormap(flipud(hot));
    xlabel('RT')
    ylabel('RE')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    xlim([min_RT max_RT]);
    ylim([min_RE max_RE]);
    pause
end

% Plot RT vs. RE per set (Rotation)
base = 15;
rot = 55;
figure;
for ii = 1:length(storage)
    RT = storage(ii).RT_set(base+1:rot);
    RE = storage(ii).RE_set(base+1:rot);
    scatter(RT,RE,30,1:length(RT),'filled');
    colormap(flipud(hot));
    xlabel('RT')
    ylabel('RE')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    xlim([min_RT max_RT]);
    ylim([min_RE max_RE]);
    pause
end

% Group average RT across trials and task
all_dm = cat(3,storage.dm_all);
all_RT = squeeze(all_dm(:,3,:));
mean_RT = nanmean(all_RT,2);
figure; plot(mean_RT); xline(120); xline(440);
xlabel('Trials'); ylabel('Mean RT')

% Group average smoothed (average per set, 8 trials per set)
set_RT = reshape(all_RT,8,[],32);
mean_RT2 = squeeze(nanmean(set_RT,1));
figure; shadedErrorBar(1:70,mean_RT2',{@nanmean,@nanstd});
xlabel('Sets'); ylabel('RT'); xline([16 56],'--r','LineWidth',2);
set(gca,'box','off','XTickLabel',[16 56],'XTick',[16 56],'YTickLabel',[],'YTick',[],...
    'FontSize',18);

% Group average RE across trials and task
all_RE = squeeze(all_dm(:,4,:));
mean_RE = nanmean(all_RE,2);
figure; plot(mean_RE); xline(120); xline(440);
xlabel('Trials'); ylabel('Mean RE');
% Group average smoothed (average per set, 8 trials per set)
set_RE = reshape(all_RE,8,[],32);
mean_RE2 = squeeze(nanmean(set_RE,1));
figure; shadedErrorBar(1:70,mean_RE2',{@nanmean,@nanstd});
xlabel(''); ylabel(''); xline([16 56],'--r','LineWidth',2);
set(gca,'box','off','XTickLabel',[16 56],'XTick',[16 56],'YTickLabel',[],'YTick',[]);


%% Balanced Integration Score (BIS)

all_task = vertcat(storage.dm_all);

% Mean, std RT across all trials and subjects
mean_RT = nanmean(all_task(:,3));
std_RT = nanstd(all_task(:,3));

% Zscore RT
for ii = 1:length(storage)
    RT = storage(ii).dm_all(:,3);
    RTz = zeros(length(RT),1);
    for jj = 1:length(RT)
        RTz(jj) = (RT(jj)-mean_RT)/std_RT;
    end
    storage(ii).RTz = RTz;
end

% Mean, std RE across all trials and subjects
mean_RE = nanmean(all_task(:,4));
std_RE = nanstd(all_task(:,4));

% Zscore RE
for ii = 1:length(storage)
    RE = storage(ii).dm_all(:,4);
    REz = zeros(length(RE),1);
    for jj = 1:length(RE)
        REz(jj) = (RE(jj)-mean_RE)/std_RE;
    end
    storage(ii).REz = REz;
end

% Add RTz and REz -> BIS for each subject
% Smaller is better
for ii = 1:length(storage)
    BIS = storage(ii).RTz + storage(ii).REz;
    storage(ii).BIS = -BIS;
end

% Plot BIS across trials for each subject
all_BIS = vertcat(storage.BIS);
min_BIS = min(all_BIS);
max_BIS = max(all_BIS);
figure;
for ii = 1:length(storage)
    scatter(1:length(storage(ii).BIS),storage(ii).BIS,20,'filled');
    xlabel("Trials")
    ylabel("BIS")
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    ylim([min_BIS max_BIS])
    xline(120)
    xline(440)
    pause
end

% Average BIS per set
all_set = storage(1).dm_all(:,5);
set_name = unique(all_set);
for ii = 1:length(storage)
    BIS = storage(ii).BIS;
    for jj = 1:length(set_name)
        BIS_set(jj) = nanmean(BIS(all_set==jj));
    end
    storage(ii).BIS_set = BIS_set;
end
set_BIS = vertcat(storage.BIS_set);
mean_BIS2 = nanmean(set_BIS,1);
figure; plot(mean_BIS2); xline(15); xline(55);
xlabel('Sets'); ylabel('Mean BIS');

% Plot BIS per set with shadedErrorBar
figure; shadedErrorBar(1:70,set_BIS,{@nanmean,@nanstd});
xlabel(''); ylabel(''); xline([16 56],'--r','LineWidth',2);
set(gca,'box','off','XTickLabel',[16 56],'XTick',[16 56],'YTickLabel',[],'YTick',[]);

% Plot BIS per set (Rotation)
all_BIS2 = horzcat(storage.BIS_set)';
min_BIS2 = min(all_BIS2,[],'all');
max_BIS2 = max(all_BIS2,[],'all');
figure;
for ii = 1:length(storage)
    BIS = storage(ii).BIS_set(16:55);
    scatter(1:length(BIS),BIS,30,'filled');
    xlabel('Set')
    ylabel('mean BIS')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    ylim([min_BIS2 max_BIS2]);
    pause
end

storage = rmfield(storage, 'RTz');
storage = rmfield(storage, 'REz');

% Check final performance of washout
valid_sub = set_BIS(:,end)>=0;
sum(valid_sub)
storage2 = storage(valid_sub==1);

% FIGURE 1C --------------------------------------------------------------
% Plot BIS per set (3rd plot)
set_BIS2 = vertcat(storage2.BIS_set);
figure; shadedErrorBar(1:70,set_BIS2,{@nanmean,@nanstd});
xlabel(''); ylabel(''); xline([15 55],'--r','LineWidth',2);
set(gca,'box','off','XTickLabel',[15 55],'XTick',[15 55],'YTickLabel',[],'YTick',[],...
    'LineWidth',1.5);

% RTz & REz set
all_rt = horzcat(storage2.RTz);
all_re = horzcat(storage2.REz);
for ii = 1:length(set_name)
    rt_set(ii,:) = nanmean(all_rt(all_set==ii,:),1);
    re_set(ii,:) = nanmean(all_re(all_set==ii,:),1);
end

% RT (1st plot)
figure; shadedErrorBar(1:70,rt_set',{@nanmean,@nanstd});
xlabel(''); ylabel(''); xline([15 55],'--r','LineWidth',2);
set(gca,'box','off','XTickLabel',[15 55],'XTick',[15 55],'YTickLabel',[],'YTick',[],...
    'LineWidth',1.5);

% RE (2nd plot)
figure; shadedErrorBar(1:70,re_set',{@nanmean,@nanstd});
xlabel(''); ylabel(''); xline([15 55],'--r','LineWidth',2);
set(gca,'box','off','XTickLabel',[15 55],'XTick',[15 55],'YTickLabel',[],'YTick',[],...
    'LineWidth',1.5);


%% Take TRs that relate to peak HRF response and are valid trials

% Create HRF
P = [6 16 1 1 6 0 32]; % Specify parameters of response function
T = 16; % Specify microtime resolution
RT = 2; % Repitition time
[hrf,~] = spm_hrf(RT,P,T); % Create hrf
% Plot HRF
figure; plot(hrf);

% Baseline + Rotation
% Create event variable for each trial
event_length = 880; % Frame length of activity
trial_num = size(trial_type_f,1); % Total number of trials
event_temp = zeros(event_length,trial_num);
counter = 1;
for ii = 1:trial_num
    event_temp(counter,ii) = 1;
    counter = counter + 2;
end
fill = zeros(8,440);
event_onset = [fill; event_temp; fill];

% Convolve each trial with hrf
for ii = 1:size(event_onset,2)
   event_hrf(:,ii) = conv(event_onset(:,ii),hrf);
end
% Plot each trial hrf
figure;
for ii = 1:size(event_hrf,2)
    plot(event_hrf(:,ii));
    hold on
    pause
end
% Match length with timeseries
event_ts = event_hrf(1:896,:);
% Find peak of each trial hrf
peak_hrf = zeros(trial_num,1);
for ii = 1:size(event_ts,2)
    [~, peak_index] = max(event_ts(:,ii));
    peak_hrf(ii) = peak_index;
end

% Find valid trials
for ii = 1:length(storage2)
    dm = storage2(ii).dm;
    valid_trial = zeros(440,1);
    valid_trial(~isnan(dm(:,3)) & ~isnan(dm(:,4))) = 1;
    dm(:,6) = valid_trial;
    storage2(ii).dm = dm;
end

% Extract valid TRs that match with peak HRF response
rot = 440;
for ii = 1:length(storage2)
    dm = storage2(ii).dm;
    BIS = storage2(ii).BIS(1:rot);
    ts = storage2(ii).ts_baserot_roi;
    ts_peak = ts(peak_hrf,:);
    valid_ts = ts_peak(dm(:,6)==1,:);
    valid_BIS = BIS(dm(:,6)==1);
    valid_dm = dm(dm(:,6)==1,:);
    storage2(ii).valid_ts = valid_ts;
    storage2(ii).valid_dm = [valid_dm valid_BIS];
end

% Washout
% Create event variable for each trial
event_length2 = 240; % Frame length of activity
trial_num2 = size(trial_type_f2,1); % Total number of trials
event_temp2 = zeros(event_length2,trial_num2);
counter = 1;
for ii = 1:trial_num2
    event_temp2(counter,ii) = 1;
    counter = counter + 2;
end
fill = zeros(8,size(event_temp2,2));
event_onset2 = [fill; event_temp2; fill];

% Convolve each trial with hrf
for ii = 1:size(event_onset2,2)
   event_hrf2(:,ii) = conv(event_onset2(:,ii),hrf);
end

% Match length with timeseries
event_ts2 = event_hrf2(1:256,:);
% Find peak of each trial hrf
peak_hrf2 = zeros(trial_num2,1);
for ii = 1:size(event_ts2,2)
    [~, peak_index2] = max(event_ts2(:,ii));
    peak_hrf2(ii) = peak_index2;
end

% Find valid trials
for ii = 1:length(storage2)
    dm = storage2(ii).dm_wash;
    valid_trial = zeros(120,1);
    valid_trial(~isnan(dm(:,3)) & ~isnan(dm(:,4))) = 1;
    dm(:,6) = valid_trial;
    storage2(ii).dm_wash = dm;
end

% Extract valid TRs that match with peak HRF response
for ii = 1:length(storage2)
    dm = storage2(ii).dm_wash;
    BIS = storage2(ii).BIS(rot+1:end);
    ts = storage2(ii).ts_wash_roi;
    ts_peak = ts(peak_hrf2,:);
    valid_ts = ts_peak(dm(:,6)==1,:);
    valid_BIS = BIS(dm(:,6)==1);
    valid_dm = dm(dm(:,6)==1,:);
    storage2(ii).valid_ts2 = valid_ts;
    storage2(ii).valid_dm2 = [valid_dm valid_BIS];
end


%% Minimum trials per set

% Baseline + Rotation
% Count how many trials across sets for each subject
sets = 1:55;
for ii = 1:length(storage2)
    dm = storage2(ii).valid_dm;
    counts = hist(dm(:,5),sets);
    set_count = [sets; counts];
    storage2(ii).set_count = set_count';
end

% Find minimum amount of trials to take from each set
counts = horzcat(storage2.set_count);
counts = counts(:,2:2:end);
min_count = min(counts,[],'all'); % Minimum trials/set = 4

% Washout
% Count how many trials across sets for each subject
sets = 56:70;
for ii = 1:length(storage2)
    dm = storage2(ii).valid_dm2;
    counts = hist(dm(:,5),sets);
    set_count = [sets; counts];
    storage2(ii).set_count2 = set_count';
end

% Find minimum amount of trials to take from each set
counts2 = horzcat(storage2.set_count2);
counts2 = counts2(:,2:2:end);
min_count2 = min(counts2,[],'all'); % Minimum trials/set = 6

min_trial = min([min_count; min_count2]);


%% Check for robustness in raw data

storage_train = storage2;

% Baseline + Rotation
% Take 4 random trials and TRs across sets for each subject
% randsample(seed, vector, amount)
for a = 1:100
    for ii = 1:length(storage_train)
        ts = storage_train(ii).valid_ts;
        dm = storage_train(ii).valid_dm;
        sets = dm(:,5);
        set_name = unique(dm(:,5));
        sample_sets = zeros(min_trial,length(set_name));
        for jj = 1:length(set_name)
            sub_trial = dm(sets==jj,1);
            sample_trials = randsample(sub_trial,min_trial,false);
            sample_sets(:,jj) = sample_trials;
        end
        [sample_sets,~] = sort(reshape(sample_sets,[],1));
        storage_train(ii).sample_sets = sample_sets;
    end
    data = horzcat(storage_train.sample_sets);
    train_sample(:,:,a) = data;
end

% Extract corresponding dm and ts data
for a = 1:size(train_sample,3)
    for ii = 1:length(storage_train)
        dm = storage_train(ii).valid_dm;
        sample_sets = train_sample(:,ii,a);
        ts = storage_train(ii).valid_ts;
        sample_dm = zeros(size(dm));
        sample_ts = zeros(size(ts));
        for jj = 1:length(sample_sets)
            sample_dm(sample_sets(jj)==dm(:,1),:) = dm(sample_sets(jj)==dm(:,1),:);
            sample_ts(sample_sets(jj)==dm(:,1),:) = ts(sample_sets(jj)==dm(:,1),:);
        end
        sample_dm(~any(sample_dm,2),:) = [];
        sample_ts(~any(sample_ts,2),:) = [];
        storage_train(ii).sample_dm = sample_dm;
        storage_train(ii).sample_ts = sample_ts;
    end
    data = vertcat(storage_train.sample_dm);
    train_dm(:,:,a) = data;
    data2 = vertcat(storage_train.sample_ts);
    train_ts(:,:,a) = data2;
end

% Washout
% Take 4 random trials and TRs across sets for each subject
% randsample(seed, vector, amount)
for a = 1:100
    for ii = 1:length(storage_train)
        ts = storage_train(ii).valid_ts2;
        dm = storage_train(ii).valid_dm2;
        sets = dm(:,5);
        set_name = unique(dm(:,5));
        sample_sets = zeros(min_trial,length(set_name));
        for jj = 1:length(set_name)
            sub_trial = dm(sets==jj+55,1); % 55 = number of sets from 1st scan
            sample_trials = randsample(sub_trial,min_trial,false);
            sample_sets(:,jj) = sample_trials;
        end
        [sample_sets,~] = sort(reshape(sample_sets,[],1));
        storage_train(ii).sample_sets2 = sample_sets;
    end
    data = horzcat(storage_train.sample_sets2);
    train_sample2(:,:,a) = data;
end

% Extract corresponding dm and ts data
for a = 1:100
    for ii = 1:length(storage_train)
        dm = storage_train(ii).valid_dm2;
        sample_sets = train_sample2(:,ii,a);
        ts = storage_train(ii).valid_ts2;
        sample_dm = zeros(size(dm));
        sample_ts = zeros(size(ts));
        for jj = 1:length(sample_sets)
            sample_dm(sample_sets(jj)==dm(:,1),:) = dm(sample_sets(jj)==dm(:,1),:);
            sample_ts(sample_sets(jj)==dm(:,1),:) = ts(sample_sets(jj)==dm(:,1),:);
        end
        sample_dm(~any(sample_dm,2),:) = [];
        sample_ts(~any(sample_ts,2),:) = [];
        storage_train(ii).sample_dm2 = sample_dm;
        storage_train(ii).sample_ts2 = sample_ts;
    end
    data = vertcat(storage_train.sample_dm2);
    data2 = vertcat(storage_train.sample_ts2);
    train_dm2(:,:,a) = data;
    train_ts2(:,:,a) = data2;
end

% sample_dm (design matrix): 
%   1 = trial number, 2 = trial type, 3 = RT, 4 = RE, 5 = set, 
%   6 = Valid, 7 = BIS
% sample_ts (timeseries): corresponding activation time-point for each
% trial
% Baseline = 60 trials (15 sets), Rotation = 160 trials (40 sets)

% Remove unnecessary fields
storage_train = rmfield(storage_train, 'dm');
storage_train = rmfield(storage_train, 'dm_wash');
storage_train = rmfield(storage_train, 'dm_all');
storage_train = rmfield(storage_train, 'ts_baserot_roi');
storage_train = rmfield(storage_train, 'ts_wash_roi');
storage_train = rmfield(storage_train, 'BIS');
storage_train = rmfield(storage_train, 'valid_ts');
storage_train = rmfield(storage_train, 'valid_ts2');
storage_train = rmfield(storage_train, 'valid_dm');
storage_train = rmfield(storage_train, 'valid_dm2');
storage_train = rmfield(storage_train, 'set_count');
storage_train = rmfield(storage_train, 'set_count2');
storage_train = rmfield(storage_train, 'sample_sets');
storage_train = rmfield(storage_train, 'sample_sets2');

% Check distribution of a single iteration
figure; histogram(train_dm(:,7,1));

% Average BIS per set across iterations
train_dmall = [train_dm; train_dm2];
train_BIS = squeeze(train_dmall(:,7,:));
set_name = train_dmall(:,5,1);
for a = 1:100
    BIS = train_BIS(:,a);
    BIS_set = zeros(70,1);
    for ii = 1:size(unique(set_name),1)
        BIS_set(ii) = mean(BIS(set_name==ii));
    end
    BIS_iter(:,a) = BIS_set;
end
figure; shadedErrorBar(1:70,mean(BIS_iter,2),std(BIS_iter,[],2));
xlabel('Sets'); ylabel('Average BIS');
xline(15); xline(55)

% Average cortical blood flow per set across iterations
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
train_tsall = [train_ts; train_ts2];
set_ids = train_dmall(:,5,1);
for a = 1:100
    ts = train_tsall(:,1:400,a);
    % Average BOLD by Yeo 7-networks
    for ii = 1:7
        mean_ts(:,ii) = mean(ts(:,schaef_id==ii),2); 
    end
    % Average BOLD by sets
    for jj = 1:length(unique(set_ids))
        meants_set(jj,:) = mean(mean_ts(set_ids==jj,:),1);
    end
    meants_iter(:,:,a) = meants_set;
end
figure; shadedErrorBar(1:70,mean(meants_iter,2),std(meants_iter,[],2));
xlabel('Sets'); ylabel('Average ts');
xline(15); xline(55);

% Cerebellar blood flow
for a = 1:100
    ts = train_tsall(:,401:428,a);
    mean_ts = mean(ts,2);
    for ii = 1:size(unique(set_name),1)
        meants_set(ii) = mean(mean_ts(set_name==ii));
    end
    meants2_iter(:,a) = meants_set;
end
figure; shadedErrorBar(1:70,mean(meants2_iter,2),std(meants2_iter,[],2));
xlabel('Sets'); ylabel('Average ts');
xline(15); xline(55);

% Basal Ganglia blood flow
for a = 1:100
    ts = train_tsall(:,429:454,a);
    mean_ts = mean(ts,2);
    for ii = 1:size(unique(set_name),1)
        meants_set(ii) = mean(mean_ts(set_name==ii));
    end
    meants3_iter(:,a) = meants_set;
end
figure; shadedErrorBar(1:70,mean(meants3_iter,2),std(meants3_iter,[],2));
xlabel('Sets'); ylabel('Average ts');
xline(15); xline(55);

% Take single iteration from random sampling for future analyses
sample_dm = train_dmall(:,:,1);
sample_ts = train_tsall(:,:,1);
base_rot = 5060;

% Reformat to align trials with subjects
trials = min_trial;
num_sets1 = 55;
num_sets2 = 15;
subjects = length(storage_train);
vars = 7;
sample1 = sample_dm(1:base_rot,:);
sample2 = sample_dm(base_rot+1:end,:);
sample1 = reshape(sample1,trials*num_sets1,subjects,vars); % trials*sets X subjects X variables
sample2 = reshape(sample2,trials*num_sets2,subjects,vars);
sample_dm = [sample1; sample2]; % trials X subject X variable

sample1 = sample_ts(1:base_rot,:);
sample2 = sample_ts(base_rot+1:end,:);
roi = 454;
sample1 = reshape(sample1,trials*num_sets1,subjects,roi);
sample2 = reshape(sample2,trials*num_sets2,subjects,roi);
sample_ts = [sample1; sample2]; % trials X subject X roi

% FIGURE 4A --------------------------------------------------------------
% Inter-trial correlation
corr_BIS2 = corr(set_BIS);
figure; imagesc(corr_BIS2);
axis('square');
set(gca,'box','off','FontSize',24,'FontName','Arial','XColor','w','YColor','w');
