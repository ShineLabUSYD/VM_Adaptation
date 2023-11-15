%% Code Outline

% The following code runs the main GLM analyses but using 1000 Schaefer cortical nodes.

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
filePattern = fullfile(myFolder, 'voltron_1000\baseline_rotation\*ses-01*.mat');
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
filePattern2 = fullfile(myFolder, 'voltron_1000\washout\*ses-01*.mat');
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
% trial type is a character array --> (n,:) gets the nth trial type

% Some subjects have trials with 0 response time
% paper removes RT>2s or RT<100ms
% Some subject have missing values for the whole set

% Outlier removal of RE
% detrending, removing trials >3SD from fitted curve (refer to paper)
% interpolation of missing trials using fitted curve

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


%% Changing outliers to NaN (alternative to removal)

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
    ts_cortex = ts(:,1:1000);
    ts_cerebellum = ts(:,1055:1082);
    ts_cnR = ts(:,1013:1020);
    ts_cnL = ts(:,1040:1047);
    ts_nacR = ts(:,1024:1025);
    ts_nacL = ts(:,1051:1052);
    ts_gpR = ts(:,1026:1027);
    ts_gpL = ts(:,1053:1054);
    ts_stn = ts(:,1099:1100);
    ts_bg = [ts_cnR ts_nacR ts_gpR ts_cnL ts_nacL ts_gpL ts_stn];
    storage(ii).ts_baserot_roi = [ts_cortex ts_cerebellum ts_bg];
end

% Washout
% Cortex, cerebellar, basal ganglia time-series (with FILL TRs)
for ii = 1:length(storage)
    ts = storage(ii).ts_wash;
    ts_cortex = ts(:,1:1000);
    ts_cerebellum = ts(:,1055:1082);
    ts_cnR = ts(:,1013:1020);
    ts_cnL = ts(:,1040:1047);
    ts_nacR = ts(:,1024:1025);
    ts_nacL = ts(:,1051:1052);
    ts_gpR = ts(:,1026:1027);
    ts_gpL = ts(:,1053:1054);
    ts_stn = ts(:,1099:1100);
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

% Group average BIS
all_BIS2 = horzcat(storage.BIS);
mean_BIS = nanmean(all_BIS2,2);
figure; plot(mean_BIS); xline(120); xline(440);
xlabel('Trials'); ylabel('Mean BIS');

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

% Average BIS per quarter (Rotation)
base = 120;
rot = 440;
set_quarter = 1:4;
set_quarter = repelem(set_quarter,80)';
for ii = 1:length(storage)
    BIS = storage(ii).BIS(base+1:rot);
    for jj = 1:length(unique(set_quarter))
        BIS_quarter(jj) = nanmean(BIS(set_quarter==jj));
    end
    storage(ii).BIS_quarter = BIS_quarter;
end
% Plot BIS per quarter
figure;
for ii = 1:length(storage)
    BIS = storage(ii).BIS_quarter;
    scatter(1:4,BIS,50,'filled');
    xlabel('Quarter')
    ylabel('mean BIS')
    plotname = storage(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    ylim([min_BIS2 max_BIS2]);
    pause
end

storage = rmfield(storage, 'RTz');
storage = rmfield(storage, 'REz');

% Check distribution of BIS for Rotation
all_BIS = horzcat(storage.BIS);
rot_BIS = all_BIS(121:440,:);
figure; histogram(rot_BIS);
% Set
all_BIS = vertcat(storage.BIS_set);
rot_BIS = all_BIS(:,16:55);
figure; histogram(rot_BIS);

% Top 50% BIS
rot_BIS = reshape(rot_BIS,[],1);
[sort_BIS,~] = sort(rot_BIS);
num = length(sort_BIS)*0.5;
top_BIS = sort_BIS(end-num+1:end);
figure; histogram(top_BIS);

% Calculate change in BIS_quarter
BIS_quarter = vertcat(storage.BIS_quarter);
delta_BIS = BIS_quarter(:,4) - BIS_quarter(:,1);
figure; histogram(delta_BIS);
figure; scatter(1:length(delta_BIS),delta_BIS,50,'filled');
yline(0); yline(1);

% Find subjects that had noticable learning or were focused at the end
% Noticable learning: delta_BIS > 1
% Focused at the end: final_BIS > 0
final_BIS = BIS_quarter(:,4);
valid_sub = delta_BIS > 1 | final_BIS > 0;
sum(valid_sub)
storage2 = storage(valid_sub==1);

% Plot BIS per quarter of new dataset
figure;
for ii = 1:length(storage2)
    BIS = storage2(ii).BIS_quarter;
    scatter(1:4,BIS,50,'filled');
    xlabel('Quarter')
    ylabel('mean BIS')
    plotname = storage2(ii).name;
    plotname = strrep(plotname,'_',' ');
    title(plotname)
    ylim([min_BIS2 max_BIS2]);
    pause
end

% Check final performance of washout
valid_sub = set_BIS(:,end)>=0;
sum(valid_sub)
storage2 = storage(valid_sub==1);

% Plot BIS per set
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
% RT
figure; shadedErrorBar(1:70,rt_set',{@nanmean,@nanstd});
xlabel(''); ylabel(''); xline([15 55],'--r','LineWidth',2);
set(gca,'box','off','XTickLabel',[15 55],'XTick',[15 55],'YTickLabel',[],'YTick',[],...
    'LineWidth',1.5);
% RE
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
% Adjust peak_hrf for different timepoints
tp = 0;
peak_hrf = peak_hrf + tp;

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
% Adjust peak hrf
tp2 = 0;
peak_hrf2 = peak_hrf2 + tp2;

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

% Split dataset in half
%storage_train = storage2(1:round(length(storage2)*0.5));
storage_train = storage2;
%storage_test = storage2(round(length(storage2)*0.5)+1:end);

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
%storage_train = rmfield(storage_train, 'BIS_quarter');
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
%base_rot = 5115; % all subjects

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
roi = 1054;
sample1 = reshape(sample1,trials*num_sets1,subjects,roi);
sample2 = reshape(sample2,trials*num_sets2,subjects,roi);
sample_ts = [sample1; sample2]; % trials X subject X roi

% Inter-trial correlation
corr_BIS2 = corr(set_BIS);
figure; imagesc(corr_BIS2);
%xline([15 55],'-w','linewidth',3);
set(gca,'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],...
    'xcolor','none','ycolor','none');



%% Peak TRs BOLD analysis

% BOLD GLM
dm_all = sample_dm;
% Swap dimensions of ts
sample_ts_rs = permute(sample_ts,[1 3 2]);

% Regress BIS against BOLD ts
coef_BIS = zeros(size(sample_ts_rs,3),size(sample_ts_rs,2)); % subject X ROI
for ii = 1:size(sample_ts_rs,3)
    sub_ts = sample_ts_rs(:,:,ii);
    glm_store_all = [];
    for jj = 1:size(sub_ts,2)
        glm = fitglm(dm_all(:,ii,7),sub_ts(:,jj));
        glm_store_all(:,jj) = glm.Coefficients.Estimate; % coef X ROI
    end
    coef_BIS(ii,:) = glm_store_all(2,:); % subject X ROI
end

% Group-level analysis (one-sample t-test)
group_coef = zeros(size(coef_BIS,2),2);
for ii = 1:size(coef_BIS,2)
    [h,p] = ttest(coef_BIS(:,ii));
    group_coef(ii,1) = h;
    group_coef(ii,2) = p;
end
roi_sig = zeros(length(group_coef),1);
mean_coef = mean(coef_BIS,1);
roi_sig(group_coef(:,1)==1) = mean_coef(group_coef(:,1)==1);
roi_pval = group_coef(group_coef(:,1)==1,2);
figure; scatter(1:length(roi_sig),roi_sig,20,'filled');

% Visualise on surface (Cortex + Cerebellum)
limits = [min(roi_sig) max(roi_sig)]; % [-0.1136 0.0978]
surf_schaef1000(roi_sig(1:1000));
surf_cbm(roi_sig(1001:1028));
% Basal Ganglia
bg_plot = zeros(1,14);
bg_plot(1,1) = mean(roi_sig(1049:1050));
bg_plot(1,3) = mean(roi_sig(1045:1048));
bg_plot(1,5) = mean(roi_sig(1051:1052));
bg_plot(1,6) = mean(roi_sig(1041:1044));
bg_plot(1,8) = mean(roi_sig(1037:1038));
bg_plot(1,10) = mean(roi_sig(1033:1036));
bg_plot(1,12) = mean(roi_sig(1039:1040));
bg_plot(1,13) = mean(roi_sig(1029:1032));
figure; plot_subcortical(bg_plot,'ventricles','False');
colormap(bluewhitered());

% Split data into 2 groups
G1 = coef_BIS(1:12,:);
G2 = coef_BIS(13:end,:);

mean_G1 = mean(G1,1);
mean_G2 = mean(G2,1);
figure; scatter(mean_G1,mean_G2,30,'filled');
[r,pval] = corr(mean_G1',mean_G2');

% One-sample t-tests
G1_coef = zeros(size(G1,2),2);
for ii = 1:size(G1,2)
    [h,p] = ttest(G1(:,ii));
    G1_coef(ii,1) = h;
    G1_coef(ii,2) = p;
end
G1_sig = zeros(length(G1_coef),1);
mean_G1 = mean(G1,1);
G1_sig(G1_coef(:,1)==1) = mean_G1(G1_coef(:,1)==1);

G2_coef = zeros(size(G2,2),2);
for ii = 1:size(G2,2)
    [h,p] = ttest(G2(:,ii));
    G2_coef(ii,1) = h;
    G2_coef(ii,2) = p;
end
G2_sig = zeros(length(G2_coef),1);
mean_G2 = mean(G2,1);
G2_sig(G2_coef(:,1)==1) = mean_G2(G2_coef(:,1)==1);

figure; scatter(G1_sig,G2_sig,30,'filled');
[r,pval] = corr(G1_sig,G2_sig);

% Index for non-significant ROIs
Ghalf_sig = G1_sig~=0 & G2_sig~=0;
Ghalf_id = zeros(size(G1_coef,1),1);
Ghalf_id(G1_coef(:,1)==1) = 1;
Ghalf_id(G2_coef(:,1)==1) = 1;

% Comparison against Neurosynth data
load('C:\PythonPrograms\NiMARE\ns_LDA400_activity.mat');
ts = squeeze(ts); % Terms X ROIs [400 X 502]
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
ns_ts = [ts_cortex ts_cerebellum ts_bg];
BISz = zscore(roi_sig); % BIS map z-scored
nsz = zscore(ns_ts'); % ROIs X Terms z-scored
ns_corr = corr(BISz,nsz);
ns_corr(isnan(ns_corr)) = 0;
% Find top 10 correlations
[~,idx] = sort(ns_corr,'descend');
% Load in neurosynth terms
[~,ns_terms,~] = xlsread('C:\PythonPrograms\NiMARE\ns_LDA400_terms.csv');
ns_terms = ns_terms(2:end,2);
ns_corr_sort = ns_corr(idx)';
ns_terms_sort = ns_terms(idx);

ns_corr_top = ns_corr_sort(1:10);
ns_terms_top = ns_terms_sort(1:10);
ns_corr_bot = ns_corr_sort(end-9:end);
ns_terms_bot = ns_terms_sort(end-9:end);
ns_corr10 = [ns_corr_top; ns_corr_bot];
ns_terms10 = vertcat(ns_terms_top,ns_terms_bot);


% Regress RT and RE against BOLD ts
coef_behav = zeros(size(sample_ts_rs,3),size(sample_ts_rs,2),2); % subject X ROI X regressor
for ii = 1:size(sample_ts_rs,3)
    sub_ts = sample_ts_rs(:,:,ii);
    glm_store_all = [];
    for jj = 1:size(sub_ts,2)
        glm = fitglm(squeeze(dm_all(:,ii,3:4)),sub_ts(:,jj));
        glm_store_all(:,jj) = glm.Coefficients.Estimate; % coef X ROI
    end
    coef_behav(ii,:,1) = glm_store_all(2,:); % subject X ROI X regressor
    coef_behav(ii,:,2) = glm_store_all(3,:);
end

% Group-level analysis (one-sample t-test)
group_coef = zeros(size(coef_behav,2),2,size(coef_behav,3));
for jj = 1:size(coef_behav,3)
    for ii = 1:size(coef_behav,2)
        [h,p] = ttest(coef_behav(:,ii,jj));
        group_coef(ii,1,jj) = h;
        group_coef(ii,2,jj) = p;
    end
end
RT_sig = zeros(size(group_coef,1),1);
RT_coef = mean(coef_behav(:,:,1),1);
RT_sig(group_coef(:,1,1)==1) = RT_coef(group_coef(:,1,1)==1);
figure; scatter(1:length(RT_sig),RT_sig,20,'filled');

RE_sig = zeros(size(group_coef,1),1);
RE_coef = mean(coef_behav(:,:,2),1);
RE_sig(group_coef(:,1,2)==1) = RE_coef(group_coef(:,1,2)==1);
figure; scatter(1:length(RE_sig),RE_sig,20,'filled');

% Visualise
% RT
limits = [min(RT_sig) max(RT_sig)];
surf_schaef1000(RT_sig(1:1000));
surf_cbm(RT_sig(401:428));
% Basal Ganglia
bg_plot = zeros(1,14);
bg_plot(1,1) = mean(RT_sig(449:450));
bg_plot(1,3) = mean(RT_sig(445:448));
bg_plot(1,5) = mean(RT_sig(451:452));
bg_plot(1,6) = mean(RT_sig(441:444));
bg_plot(1,8) = mean(RT_sig(437:438));
bg_plot(1,10) = mean(RT_sig(433:436));
bg_plot(1,12) = mean(RT_sig(439:440));
bg_plot(1,13) = mean(RT_sig(429:432));
figure; plot_subcortical(bg_plot,'ventricles','False');
colormap(bluewhitered());

% RE
limits = [min(RE_sig) max(RE_sig)]; %[-0.0093 0.0051]
surf_schaef1000(RE_sig(1:1000));
surf_cbm(RE_sig(1001:1028));
% Basal Ganglia
bg_plot = zeros(1,14);
bg_plot(1,1) = mean(RE_sig(449:450));
bg_plot(1,3) = mean(RE_sig(445:448));
bg_plot(1,5) = mean(RE_sig(451:452));
bg_plot(1,6) = mean(RE_sig(441:444));
bg_plot(1,8) = mean(RE_sig(437:438));
bg_plot(1,10) = mean(RE_sig(433:436));
bg_plot(1,12) = mean(RE_sig(439:440));
bg_plot(1,13) = mean(RE_sig(429:432));
figure; plot_subcortical(bg_plot,'ventricles','False');
colormap(bluewhitered());

[r,pval] = corr(mean_coef',-RT_coef');
[r,pval] = corr(mean_coef',-RE_coef');

% Regress Task against BOLD ts

% Create task regressor matrix
% Create task regressor matrix
trial_label = dm_all(:,1,2);
task_dm = zeros(280,3);
task_dm(61:140,1) = 1;
task_dm(141:220,2) = 1;
task_dm(trial_label==2,3) = 1;

% Fit GLM with task_dm
coef_task = zeros(size(sample_ts_rs,3),size(sample_ts_rs,2),4); % subject X ROI X regressor
for ii = 1:size(sample_ts_rs,3)
    sub_ts = sample_ts_rs(:,:,ii);
    glm_store_all = [];
    for jj = 1:size(sub_ts,2)
        glm = fitglm(task_dm,sub_ts(:,jj));
        glm_store_all(:,jj) = glm.Coefficients.Estimate; % coef X ROI
    end
    coef_task(ii,:,1) = glm_store_all(1,:); % base
    coef_task(ii,:,2) = glm_store_all(2,:); % rot early
    coef_task(ii,:,3) = glm_store_all(3,:); % rot late
    coef_task(ii,:,4) = glm_store_all(4,:); % washout
end

% Group one-sample t-test
group_coef = zeros(size(coef_task,2),2,size(coef_task,3));
for jj = 1:size(coef_task,3)
    for ii = 1:size(coef_task,2)
        [h,p] = ttest(coef_task(:,ii,jj));
        group_coef(ii,1,jj) = h;
        group_coef(ii,2,jj) = p; % ROI X stat X task
    end
end

% Base ROIs
base_sig = zeros(size(group_coef,1),1);
mean_base = mean(coef_task(:,:,1),1);
base_sig(group_coef(:,1,1)==1) = mean_base(group_coef(:,1,1)==1);
figure; scatter(1:length(base_sig),base_sig,20,'filled');

% Visualise on surface (Cortex + Cerebellum)
surf_schaef2(base_sig(1:400));
surf_cbm(base_sig(401:428));
% Basal Ganglia
bg_plot = zeros(1,14);
bg_plot(1,1) = mean(base_sig(449:450));
bg_plot(1,3) = mean(base_sig(445:448));
bg_plot(1,5) = mean(base_sig(451:452));
bg_plot(1,6) = mean(base_sig(441:444));
bg_plot(1,8) = mean(base_sig(437:438));
bg_plot(1,10) = mean(base_sig(433:436));
bg_plot(1,12) = mean(base_sig(439:440));
bg_plot(1,13) = mean(base_sig(429:432));
figure; plot_subcortical(bg_plot,'ventricles','False');
colormap(bluewhitered());


% Rotation ROIs
rotE_sig = zeros(size(group_coef,1),1);
mean_rotE = mean(coef_task(:,:,2),1);
rotE_sig(group_coef(:,1,2)==1) = mean_rotE(group_coef(:,1,2)==1);
figure; scatter(1:length(rotE_sig),rotE_sig,20,'filled');

% Visualise on surface (Cortex + Cerebellum)
surf_schaef2(rotE_sig(1:400));
surf_cbm(rotE_sig(401:428));
% Basal Ganglia
bg_plot = zeros(1,14);
bg_plot(1,1) = mean(rotE_sig(449:450));
bg_plot(1,3) = mean(rotE_sig(445:448));
bg_plot(1,5) = mean(rotE_sig(451:452));
bg_plot(1,6) = mean(rotE_sig(441:444));
bg_plot(1,8) = mean(rotE_sig(437:438));
bg_plot(1,10) = mean(rotE_sig(433:436));
bg_plot(1,12) = mean(rotE_sig(439:440));
bg_plot(1,13) = mean(rotE_sig(429:432));
figure; plot_subcortical(bg_plot,'ventricles','False');
colormap(bluewhitered());

% Rot late
rotL_sig = zeros(size(group_coef,1),1);
mean_rotL = mean(coef_task(:,:,3),1);
rotL_sig(group_coef(:,1,3)==1) = mean_rotL(group_coef(:,1,3)==1);
figure; scatter(1:length(rotL_sig),rotL_sig,20,'filled');

% Visualise on surface (Cortex + Cerebellum)
surf_schaef2(rotL_sig(1:400));
surf_cbm(rotL_sig(401:428));
% Basal Ganglia
bg_plot = zeros(1,14);
bg_plot(1,1) = mean(rotL_sig(449:450));
bg_plot(1,3) = mean(rotL_sig(445:448));
bg_plot(1,5) = mean(rotL_sig(451:452));
bg_plot(1,6) = mean(rotL_sig(441:444));
bg_plot(1,8) = mean(rotL_sig(437:438));
bg_plot(1,10) = mean(rotL_sig(433:436));
bg_plot(1,12) = mean(rotL_sig(439:440));
bg_plot(1,13) = mean(rotL_sig(429:432));
figure; plot_subcortical(bg_plot,'ventricles','False');
colormap(bluewhitered());

% Washout ROIs
wash_sig = zeros(size(group_coef,1),1);
mean_wash = mean(coef_task(:,:,4),1);
wash_sig(group_coef(:,1,4)==1) = mean_wash(group_coef(:,1,4)==1);
figure; scatter(1:length(wash_sig),wash_sig,20,'filled');

% Visualise on surface (Cortex + Cerebellum)
surf_schaef2(wash_sig(1:400));
surf_cbm(wash_sig(401:428));
% Basal Ganglia
bg_plot = zeros(1,14);
bg_plot(1,1) = mean(wash_sig(449:450));
bg_plot(1,3) = mean(wash_sig(445:448));
bg_plot(1,5) = mean(wash_sig(451:452));
bg_plot(1,6) = mean(wash_sig(441:444));
bg_plot(1,8) = mean(wash_sig(437:438));
bg_plot(1,10) = mean(wash_sig(433:436));
bg_plot(1,12) = mean(wash_sig(439:440));
bg_plot(1,13) = mean(wash_sig(429:432));
figure; plot_subcortical(bg_plot,'ventricles','False');
colormap(bluewhitered());

% Compare task maps
figure;
scatter(base_sig,rot_sig,30,'filled');

figure;
scatter(base_sig,wash_sig,30,'filled');

figure;
scatter(wash_sig,rot_sig,30,'filled');

% PCA of task betas
coef_task_rs = permute(coef_task,[1 3 2]); % subject X task X ROI
coef_task_rs = reshape(coef_task_rs,[],454);
[task_vec,task_val,~,~,task_explained,~] = pca(coef_task_rs);
figure; scatter(1:length(task_explained),task_explained,30,'filled');
yline(1);

% Regress BIS by task (4 tasks)
% Create task labels
trial_label = dm_all(:,1,2);
task_id = zeros(280,1);
task_id(trial_label==0) = 0;
task_id(61:140) = 1;
task_id(141:220) = 2;
task_id(trial_label==2) = 3;

% Baseline
base_BIS = zeros(size(sample_ts_rs,3),size(sample_ts_rs,2)); % subject X ROI
for ii = 1:size(sample_ts_rs,3)
    sub_ts = sample_ts_rs(task_id==0,:,ii);
    glm_store_all = [];
    for jj = 1:size(sub_ts,2)
        glm = fitglm(dm_all(task_id==0,ii,7),sub_ts(:,jj));
        glm_store_all(:,jj) = glm.Coefficients.Estimate; % coef X ROI
    end
    base_BIS(ii,:) = glm_store_all(2,:); % subject X ROI
end

% Group-level analysis (one-sample t-test)
group_coef = zeros(size(base_BIS,2),2);
for ii = 1:size(base_BIS,2)
    [h,p] = ttest(base_BIS(:,ii));
    group_coef(ii,1) = h;
    group_coef(ii,2) = p;
end
base_sig2 = zeros(length(group_coef),1);
mean_coef = mean(base_BIS,1);
base_sig2(group_coef(:,1)==1) = mean_coef(group_coef(:,1)==1);
figure; scatter(1:length(base_sig2),base_sig2,20,'filled');

% Early Rotation
rotE_BIS = zeros(size(sample_ts_rs,3),size(sample_ts_rs,2)); % subject X ROI
for ii = 1:size(sample_ts_rs,3)
    sub_ts = sample_ts_rs(task_id==1,:,ii);
    glm_store_all = [];
    for jj = 1:size(sub_ts,2)
        glm = fitglm(dm_all(task_id==1,ii,7),sub_ts(:,jj));
        glm_store_all(:,jj) = glm.Coefficients.Estimate; % coef X ROI
    end
    rotE_BIS(ii,:) = glm_store_all(2,:); % subject X ROI
end

% Group-level analysis (one-sample t-test)
group_coef = zeros(size(rotE_BIS,2),2);
for ii = 1:size(rotE_BIS,2)
    [h,p] = ttest(rotE_BIS(:,ii));
    group_coef(ii,1) = h;
    group_coef(ii,2) = p;
end
rotE_sig2 = zeros(length(group_coef),1);
mean_coef = mean(rotE_BIS,1);
rotE_sig2(group_coef(:,1)==1) = mean_coef(group_coef(:,1)==1);
figure; scatter(1:length(rotE_sig2),rotE_sig2,20,'filled');

% Late Rotation
rotL_BIS = zeros(size(sample_ts_rs,3),size(sample_ts_rs,2)); % subject X ROI
for ii = 1:size(sample_ts_rs,3)
    sub_ts = sample_ts_rs(task_id==2,:,ii);
    glm_store_all = [];
    for jj = 1:size(sub_ts,2)
        glm = fitglm(dm_all(task_id==2,ii,7),sub_ts(:,jj));
        glm_store_all(:,jj) = glm.Coefficients.Estimate; % coef X ROI
    end
    rotL_BIS(ii,:) = glm_store_all(2,:); % subject X ROI
end

% Group-level analysis (one-sample t-test)
group_coef = zeros(size(rotL_BIS,2),2);
for ii = 1:size(rotL_BIS,2)
    [h,p] = ttest(rotL_BIS(:,ii));
    group_coef(ii,1) = h;
    group_coef(ii,2) = p;
end
rotL_sig2 = zeros(length(group_coef),1);
mean_coef = mean(rotL_BIS,1);
rotL_sig2(group_coef(:,1)==1) = mean_coef(group_coef(:,1)==1);
figure; scatter(1:length(rotL_sig2),rotL_sig2,20,'filled');

% Washout
wash_BIS = zeros(size(sample_ts_rs,3),size(sample_ts_rs,2)); % subject X ROI
for ii = 1:size(sample_ts_rs,3)
    sub_ts = sample_ts_rs(task_id==3,:,ii);
    glm_store_all = [];
    for jj = 1:size(sub_ts,2)
        glm = fitglm(dm_all(task_id==3,ii,7),sub_ts(:,jj));
        glm_store_all(:,jj) = glm.Coefficients.Estimate; % coef X ROI
    end
    wash_BIS(ii,:) = glm_store_all(2,:); % subject X ROI
end

% Group-level analysis (one-sample t-test)
group_coef = zeros(size(wash_BIS,2),2);
for ii = 1:size(wash_BIS,2)
    [h,p] = ttest(wash_BIS(:,ii));
    group_coef(ii,1) = h;
    group_coef(ii,2) = p;
end
wash_sig2 = zeros(length(group_coef),1);
mean_coef = mean(wash_BIS,1);
wash_sig2(group_coef(:,1)==1) = mean_coef(group_coef(:,1)==1);
figure; scatter(1:length(wash_sig2),wash_sig2,20,'filled');

% Compare task BIS maps
mean_base = mean(base_BIS,1);
mean_rotE = mean(rotE_BIS,1);
mean_rotL = mean(rotL_BIS,1);
mean_wash = mean(wash_BIS,1);

task_mean = vertcat(mean_base,mean_rotE,mean_rotL,mean_wash)';
[task_corr,task_pval] = corr(task_mean);
figure; imagesc(task_corr);

task_lower = tril(task_corr,-1);

% Spin test correlations
schaef_centroids = readtable(['C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\' ...
    'Functions\schaefer_parcellation\' ...
    'Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'],'Range','C2:E1001');
schaef_centroids = table2array(schaef_centroids);
perm_id1000 = rotate_parcellation(schaef_centroids(1:500,:),schaef_centroids(501:end,:),5000);
%load('C:\Users\JoshB\OneDrive\Documents\University\PhD\visuomotor_Gallivan\Code\data\perm_id1000.mat');
p_perm = perm_sphere_p(mean_G1(1:1000)',mean_G2(1:1000)',perm_id1000,"pearson");
