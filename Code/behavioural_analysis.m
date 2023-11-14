%% Code Outline

% This code produces the results for behavioural analyses


%% Statistical Comparisons of Behavioural data (BIS)

% Load in 'sample_dm.mat' from 'dataprep.mat'

all_BIS = sample_dm(:,:,7);
set_id = sample_dm(:,:,5);

% BASELINE - Early vs. Late
base_first = all_BIS(set_id(:,1)==1,:);
base_last = all_BIS(set_id(:,1)==15,:);

mean_bf = nanmean(base_first,1);
mean_bl = nanmean(base_last,1);

% Boxplots
figure; boxplot([mean_bf' mean_bl'],'Labels',{'Set 1','Set 15'},'Colors','k');
ylabel('Mean BIS')
hold on
% Joined scatter plot
N = length(storage2);
group_BIS = [mean_bf'; mean_bl'];
x =[repelem(1,N); repelem(2,N)];
x = x';
x = reshape(x,2*N,1);
%figure
scatter(x,group_BIS,20,'black','filled','jitter','on','jitterAmount',0.05);
line([x(1:N) x(N+1:end)]',[group_BIS(1:N) group_BIS(N+1:end)]','Color','black')
xlim([0.5 2.5])

% Paired sample t-test
[h,p,ci,stats] = ttest(mean_bf',mean_bl');
mean_diff = mean(mean_bl) - mean(mean_bf);
mean_std = std([mean_bf';mean_bl']);


% ROTATION - Early vs. Late
rot_first = all_BIS(set_id(:,1)==16,:);
rot_last = all_BIS(set_id(:,1)==55,:);

mean_rf = nanmean(rot_first,1);
mean_rl = nanmean(rot_last,1);

% Boxplots
figure; boxplot([mean_rf' mean_rl'],'Labels',{'Set 16','Set 55'});
ylabel('Mean BIS')
hold on
% Joined scatter plot
N = length(storage2);
group_BIS = [mean_rf'; mean_rl'];
x =[repelem(1,N); repelem(2,N)];
x = x';
x = reshape(x,2*N,1);
%figure
scatter(x,group_BIS)
line([x(1:N) x(N+1:end)]',[group_BIS(1:N) group_BIS(N+1:end)]')
xlim([0.5 2.5])

% Paired sample t-test
[h,p,ci,stats] = ttest(mean_rf',mean_rl');
mean_diff = mean(mean_rl) - mean(mean_rf);
mean_std = std([mean_rf';mean_rl']);


% WASHOUT - Early vs. Late
wash_first = all_BIS(set_id(:,1)==56,:);
wash_last = all_BIS(set_id(:,1)==70,:);

mean_wf = nanmean(wash_first,1);
mean_wl = nanmean(wash_last,1);

% Boxplots
figure; boxplot([mean_wf' mean_wl'],'Labels',{'Set 56','Set 70'});
ylabel('Mean BIS')
hold on
% Joined scatter plot
N = length(storage2);
group_BIS = [mean_wf'; mean_wl'];
x =[repelem(1,N); repelem(2,N)];
x = x';
x = reshape(x,2*N,1);
%figure
scatter(x,group_BIS)
line([x(1:N) x(N+1:end)]',[group_BIS(1:N) group_BIS(N+1:end)]')
xlim([0.5 2.5])

% Paired sample t-test
[h,p,ci,stats] = ttest(mean_wf',mean_wl');
mean_diff = mean(mean_wl) - mean(mean_wf);
mean_std = std([mean_wf';mean_wl']);


% Late Baseline vs. Late Rotation vs. Late Washout (peak performance)
% FIGURE 2B --------------------------------------------------------------
% Box plot with 3 groups
data1 = mean_bl';
data2 = mean_rl';
data3 = mean_wl';
% Colormap
RGB_color = [0 63 92;188 80 144;255 166 0]/255;
RGB_color2(:,1) = repelem(RGB_color(:,1),23);
RGB_color2(:,2) = repelem(RGB_color(:,2),23);
RGB_color2(:,3) = repelem(RGB_color(:,3),23);
% Group labels
N = length(data1);
group_BIS = [data1; data2; data3];
x =[repelem(1,N); repelem(2,N); repelem(3,N)];
x = x';
x = reshape(x,3*N,1);
% Box plot
figure; boxplot([data1 data2 data3],'Labels',{'','',''},'Colors',RGB_color,'Symbol','','ColorGroup',x);
ylim([-1.5 3]);
hold on
% Scatter plot
scatter(x,group_BIS,20,RGB_color2,'filled','jitter','on','jitterAmount',0.1);
set(gca,'FontSize',24,'FontName','Arial','linew',1.5);
set(findobj(gca,'type','line'),'linew',2);

% One-way Anova
task_id = {'Baseline','Rotation','Washout'}';
task = repelem(task_id,length(mean_bl));
late_matrix = [mean_bl'; mean_rl'; mean_wl'];
[p,tbl,stats] = anova1(late_matrix,task);


% Early Baseline vs. Early Rotation vs. Early Washout (starting
% performance)
% One-way Anova
task_id = {'Baseline','Rotation','Washout'}';
task = repelem(task_id,22);
early_matrix = [mean_bf'; mean_rf'; mean_wf'];
[p,tbl,stats] = anova1(early_matrix,task);

% Paired sample t-test
[h,p,ci,stats] = ttest(mean_bf',mean_wf');


% Baseline learning vs. Rotation learning vs. Washout learning
base_learn = mean_bl - mean_bf;
rot_learn = mean_rl - mean_rf;
wash_learn = mean_wl - mean_wf;
% One-way Anova
learn_matrix = [base_learn'; rot_learn'; wash_learn'];
[p,tbl,stats] = anova1(learn_matrix,task);

% FIGURE 2C --------------------------------------------------------------
% Box chart
data = [mean_bf'; mean_bl'; mean_rf'; mean_rl'; mean_wf'; mean_wl'];
group = [1; 1; 3; 3; 5; 5];
group = repelem(group,23);
cgroup = [1; 2; 1; 2; 1; 2];
cgroup = repelem(cgroup,23,1);

figure;
b = boxchart(group,data,'GroupByColor',cgroup,'MarkerStyle','none');
b(1).BoxFaceColor = [0 63 92]./255;
b(2).BoxFaceColor = [255 166 0]./255;

hold on
RGB_color = [0 63 92; 255 166 0]./255;
RGB_color2(:,1) = repelem(RGB_color(:,1),23);
RGB_color2(:,2) = repelem(RGB_color(:,2),23);
RGB_color2(:,3) = repelem(RGB_color(:,3),23);
RGB_color3 = repmat(RGB_color2,3,1);
group2 = [0.75; 1.25; 2.75; 3.25; 4.75; 5.25];
group2 = repelem(group2,23);
scatter(group2,data,20,RGB_color3,'filled','jitter','on','jitterAmount',0.1);

ylim([-4 4]);
set(gca,'FontSize',24,'FontName','Arial','linew',1.5, ...
    'XTick',[1 3 5],'XTickLabel',[],'box','off');
set(findobj(gca,'type','line'),'linew',2);