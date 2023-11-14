%% Code Outline

% surf_schaef2 = cerebral cortex plot
% surf_cbm = cerebellum plot
% subcort_plot or plot_subcortical = subcortical plots


%% Peak TRs BOLD analysis

% Load in 'sample_dm.mat' and 'sample_ts.mat' from dataprep

% Check mean BOLD maps
mean_ts = squeeze(mean(sample_ts,[1 2]));
limits = [min(mean_ts) max(mean_ts)]; % [-0.0606, 0.0472]
surf_schaef2(mean_ts(1:400));
surf_cbm(mean_ts(401:428));
subcort_plot(mean_ts);colormap(multigradient('preset','div.cb.spectral.10'));

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
limits = [min(roi_sig) max(roi_sig)]; % [-0.1118 0.0749];
surf_schaef2(roi_sig(1:400));
surf_cbm(roi_sig(401:428));
% Basal Ganglia
subcort_plot(roi_sig);
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

% Spin Permutations
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\statistic_testing\perm_id.mat');
[p_perm, null_dist] = perm_sphere_p(mean_coef(1:400)',RT_coef(1:400)',perm_id,"pearson");

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
RT_sig = -1.*RT_sig;
limits = [min(RT_sig) max(RT_sig)]; % [-0.7414 0.1696];
surf_schaef2(RT_sig(1:400));
surf_cbm(RT_sig(401:428));
% Basal Ganglia
subcort_plot(RT_sig);
colormap(bluewhitered());

% RE
RE_sig = -1.*RE_sig;
limits = [min(RE_sig) max(RE_sig)]; % [-0.0051 0.0069];
surf_schaef2(RE_sig(1:400));
surf_cbm(RE_sig(401:428));
% Basal Ganglia
subcort_plot(RE_sig);
colormap(bluewhitered());

% Compare RT, RE and BIS betas
[rho,pval] = corr(mean_coef',-RT_coef');
[rho,pval] = corr(mean_coef',-RE_coef');

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
task_lower = tril(task_corr,-1);


%% Figures

% FIGURE 3C --------------------------------------------------------------
data1 = mean_coef;
data2 = -RT_coef;
cmap = 'black';
figure; 
scatter(data1,data2,50,cmap,'filled');
h1 = lsline();
h1.Color = 'r';
h1.LineWidth = 2;
set(gca,'box','off','FontSize',24,'FontName','Arial','linew',1.5);

data1 = mean_coef;
data2 = -RE_coef;
cmap = 'black';
figure; 
scatter(data1,data2,50,cmap,'filled');
h1 = lsline();
h1.Color = 'r';
h1.LineWidth = 2;
set(gca,'box','off','FontSize',24,'FontName','Arial','linew',1.5);


% FIGURE 3D --------------------------------------------------------------
% Heatmap
data = task_lower;
data = round(task_lower,2);
data(data==0) = nan;

figure;
heatmap(data); 
set(gca,'MissingDataColor','w','FontSize',18,'FontName','Arial');


% FIGURE 3E --------------------------------------------------------------
% Scatter plot with transparency
group_id = Ghalf_sig;

data1 = mean_G1;
data2 = mean_G2;

data_sig = [mean_G1(Ghalf_sig==1)' mean_G2(Ghalf_sig==1)'];
data_nsig = [mean_G1(Ghalf_sig==0)' mean_G2(Ghalf_sig==0)'];

cmap = zeros(size(data_nsig,1),3);
cmap(:,1) = 175/255;
cmap(:,2) = 171/255;
cmap(:,3) = 171/255;

figure; 
s2 = scatter(data_nsig(:,1),data_nsig(:,2),50,cmap,'filled');
s2.MarkerFaceAlpha = 0.3;
hold on
s1 = scatter(data_sig(:,1),data_sig(:,2),50,'black','filled');

h = lsline();
h(1).Color = 'r';
h(1).LineWidth = 2;
h(2).Color = 'none';
set(gca,'box','off','FontSize',24,'FontName','Arial','linew',1.5);
hold off


% FIGURE 3F --------------------------------------------------------------
% Top panel
data = ns_corr_top;
data = data(1:5)';

min_val = 0;
max_val = 0.6;
figure;
spider_plot(data, ...
    'AxesLimits',[repelem(min_val,5); repelem(max_val,5)], ...
    'AxesInterval',3, ...
    'AxesDisplay','one',...
    'FillOption','on', ...
    'FillTransparency',0.3, ...
    'AxesRadial','on',...
    'AxesFont','Arial', ...
    'LabelFont','Arial', ...
    'AxesFontSize',16, ...
    'LineWidth',2, ...
    'MinorGrid','off', ...
    'Color',[1 0 0]);

% Bottomg panel
data = abs(flipud(ns_corr_bot)');
data = data(1:5)';

min_val = 0;
max_val = 0.6;
figure;
spider_plot(data, ...
    'AxesLimits',[repelem(min_val,5); repelem(max_val,5)], ...
    'AxesInterval',3, ...
    'AxesDisplay','one',...
    'FillOption','on', ...
    'FillTransparency',0.3, ...
    'AxesRadial','on',...
    'AxesFont','Arial', ...
    'LabelFont','Arial', ...
    'AxesFontSize',16, ...
    'LineWidth',2, ...
    'MinorGrid','off', ...
    'Color',[0 0 1]);