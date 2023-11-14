%% Code Outline

% The following code:
%       - Applies k-means clustering to BIS scores
%       - Calculates the correlation scores for each cluster and task
%       condition
%       - Calculates energy landscape for each cluster and task
%       condition


%% Group clustering

% Load in 'sample_dm.mat' and 'sample_ts.mat' from 'dataprep.mat'

% Clustering groups based off performance
min_trial = 4;
all_BIS = sample_dm(:,:,7);
num_sets = 70;
trials = min_trial;
subjects = size(sample_dm,2);
all_BIS2 = reshape(all_BIS,trials,num_sets,subjects);
set_BIS = squeeze(nanmean(all_BIS2,1))';
% K-means clustering of set_BIS
k = 2; % Number of clusters
BIS_kmeans = []; % Store clustering results
for ii = 1:subjects-2
    for jj = 1:100
        [idx,C] = kmeans(set_BIS,k,'MaxIter',100);
        BIS_kmeans(:,jj,ii) = idx; % subject X iteration X cluster size
    end
    k = k + 1;
end

% Adjusted Mutual Information (AMI)
% Run AMI for each iteration against every other iteration within a k value
clusters = size(BIS_kmeans,3);
BIS_ami = zeros(100,100,clusters);
for ii = 1:size(BIS_kmeans,3)
    sim_all = zeros(100,100);
    for jj = 1:size(BIS_kmeans,2)
       for kk = 1:size(BIS_kmeans,2)
           temp = ami(BIS_kmeans(:,jj,ii),BIS_kmeans(:,kk,ii));
           sim_all(kk,jj) = temp;
       end  
    end
    BIS_ami(:,:,ii) = sim_all; % iteration X iteration X cluster size
end

% FIGURE 4B --------------------------------------------------------------
% Plot mean and SD of AMI for cluster stability
BIS_ami2 = reshape(BIS_ami,100^2,clusters);
figure; shadedErrorBar(2:1:size(BIS_kmeans,1),mean(BIS_ami2),std(BIS_ami2,[],1));
xlabel('Cluster sizes'); ylabel('Mean AMI');

% Rerun kmeans clustering with most stable cluster size
% Number of clusters = 3
k = 3;
[BIS_idx,C] = kmeans(set_BIS,k,'MaxIter',100);
% Exploration of clusters
BIS_c1 = set_BIS(BIS_idx==1,:);
BIS_c2 = set_BIS(BIS_idx==2,:);
BIS_c3 = set_BIS(BIS_idx==3,:);

min_BIS = min(set_BIS,[],'all');
max_BIS = max(set_BIS,[],'all');

% Individual plots
subplot(1,3,1)
plot(mean(BIS_c1));
xlabel('Sets'); ylabel('BIS'); title('C1');
xline([16 56],'LineWidth',2); 
axis tight
ylim([min_BIS max_BIS]);

subplot(1,3,2)
plot(mean(BIS_c2,1));
xlabel('Sets'); ylabel('BIS'); title('C2');
xline([16 56],'LineWidth',2); 
axis tight
ylim([min_BIS max_BIS]);

subplot(1,3,3)
plot(mean(BIS_c3,1));
xlabel('Sets'); ylabel('BIS'); title('C3');
xline([16 56],'LineWidth',2); 
axis tight
ylim([min_BIS max_BIS]);


% FIGURE 4C --------------------------------------------------------------
% Multiple line plots with shaded error bar (3 groups)
data1 = BIS_c1;
data2 = BIS_c2;
data3 = BIS_c3;

figure;
shadedErrorBar(1:70,mean(data1,1),std(data1,[],1)/sqrt(size(data1,1)),{'-','color',[222 177 97]./255,'LineWidth',1.5},0);

hold on
shadedErrorBar(1:70,mean(data2,1),std(data2,[],1)/sqrt(size(data2,1)),{'-','color',[53 155 124]./255,'LineWidth',1.5},0.3);
shadedErrorBar(1:70,mean(data3,1),std(data3,[],1)/sqrt(size(data3,1)),{'-','color',[34 102 141]./255,'LineWidth',1.5},0.3);

xline([16 56],'--r','LineWidth',2);
%set(gca,'box','off','XTickLabel',[16 56],'XTick',[16 56],'YTickLabel',[0],'YTick',[0],...
    %'FontSize',24, 'FontName', 'Arial','linew',1.5);
set(gca,'box','off','XTickLabel',[16 56],'XTick',[16 56],'FontSize',24,'FontName','Arial','linew',1.5);
axis('tight');



% Mean BIS across sets
mean_C1 = mean(BIS_c1,1);
mean_C2 = mean(BIS_c2,1);
mean_C3 = mean(BIS_c3,1);

% Time-series for each group
ts_c1 = sample_ts(:,BIS_idx==1,:);
ts_c2 = sample_ts(:,BIS_idx==2,:);
ts_c3 = sample_ts(:,BIS_idx==3,:);

% Average baseline BOLD
task_label = squeeze(sample_dm(:,1,2));
base_c1 = squeeze(mean(ts_c1(task_label==0,:,:),[1 2]));
base_c2 = squeeze(mean(ts_c2(task_label==0,:,:),[1 2]));
base_c3 = squeeze(mean(ts_c3(task_label==0,:,:),[1 2]));

surf_schaef2(base_c1(1:400));
surf_cbm(base_c1(401:428));

surf_schaef2(base_c2(1:400));
surf_cbm(base_c2(401:428));

surf_schaef2(base_c3(1:400));
surf_cbm(base_c3(401:428));

base_diff = base_c2 - base_c3;
surf_schaef2(base_diff(1:400));
surf_cbm(base_diff(401:428));
subcort_plot(base_diff); colormap(bluewhitered());

rot_c1 = squeeze(mean(ts_c1(task_label==1,:,:),[1 2]));
rot_c2 = squeeze(mean(ts_c2(task_label==1,:,:),[1 2]));
rot_c3 = squeeze(mean(ts_c3(task_label==1,:,:),[1 2]));

rot_diff = rot_c2 - rot_c3;
surf_schaef2(rot_diff(1:400));
surf_cbm(rot_diff(401:428));
subcort_plot(rot_diff); colormap(bluewhitered());

wash_c1 = squeeze(mean(ts_c1(task_label==2,:,:),[1 2]));
wash_c2 = squeeze(mean(ts_c2(task_label==2,:,:),[1 2]));
wash_c3 = squeeze(mean(ts_c3(task_label==2,:,:),[1 2]));

wash_diff = wash_c2 - wash_c3;
surf_schaef2(wash_diff(1:400));
surf_cbm(wash_diff(401:428));
subcort_plot(wash_diff); colormap(bluewhitered());

surf_schaef2(mean_diff(1:400));
surf_cbm(mean_diff(401:428));
subcort_plot(mean_diff); colormap(bluewhitered());

% Mean BOLD for each cluster and correlate BOLD across trials
c1_mean = squeeze(mean(ts_c1,2));
c2_mean = squeeze(mean(ts_c2,2));
c3_mean = squeeze(mean(ts_c3,2));

% Separate mean ts by task
c1_base = c1_mean(1:60,:);
c2_base = c2_mean(1:60,:);
c3_base = c3_mean(1:60,:);

c1_rotE = c1_mean(61:140,:);
c2_rotE = c2_mean(61:140,:);
c3_rotE = c3_mean(61:140,:);

c1_rotL = c1_mean(141:220,:);
c2_rotL = c2_mean(141:220,:);
c3_rotL = c3_mean(141:220,:);

c1_wash = c1_mean(221:280,:);
c2_wash = c2_mean(221:280,:);
c3_wash = c3_mean(221:280,:);


% Calculate correlation between brain maps
c1_base_corr = corr(c1_base');
c2_base_corr = corr(c2_base');
c3_base_corr = corr(c3_base');

c1_rotE_corr = corr(c1_rotE');
c2_rotE_corr = corr(c2_rotE');
c3_rotE_corr = corr(c3_rotE');

c1_rotL_corr = corr(c1_rotL');
c2_rotL_corr = corr(c2_rotL');
c3_rotL_corr = corr(c3_rotL');

c1_wash_corr = corr(c1_wash');
c2_wash_corr = corr(c2_wash');
c3_wash_corr = corr(c3_wash');

% Get lower triangle
base_temp = tril(ones(60)-eye(60));
base_id = find(base_temp);

rot_temp = tril(ones(80)-eye(80));
rot_id = find(rot_temp);

c1_base_low = c1_base_corr(base_id);
c2_base_low = c2_base_corr(base_id);
c3_base_low = c3_base_corr(base_id);
c1_wash_low = c1_wash_corr(base_id);
c2_wash_low = c2_wash_corr(base_id);
c3_wash_low = c3_wash_corr(base_id);

c1_rotE_low = c1_rotE_corr(rot_id);
c2_rotE_low = c2_rotE_corr(rot_id);
c3_rotE_low = c3_rotE_corr(rot_id);
c1_rotL_low = c1_rotL_corr(rot_id);
c2_rotL_low = c2_rotL_corr(rot_id);
c3_rotL_low = c3_rotL_corr(rot_id);

% Explore
figure; histogram(abs(c1_base_low));
hold on
histogram(abs(c2_base_low));
histogram(abs(c3_base_low));

figure; histogram(abs(c1_rotE_low));
hold on
histogram(abs(c2_rotE_low));
histogram(abs(c3_rotE_low));

figure; histogram(abs(c1_rotL_low));
hold on
histogram(abs(c2_rotL_low));
histogram(abs(c3_rotL_low));

figure; histogram(abs(c1_wash_low));
hold on
histogram(abs(c2_wash_low));
histogram(abs(c3_wash_low));

% Statistic testing
[h,p,D] = kstest2(c1_rotE_low,c2_rotE_low);
[h,p,D] = kstest2(c1_rotE_low,c3_rotE_low);
[h,p,D] = kstest2(c2_rotE_low,c3_rotE_low);

[h,p,D] = kstest2(c2_rotL_low,c3_rotL_low);
[h,p,D] = kstest2(c1_rotL_low,c3_rotL_low);
[h,p,D] = kstest2(c1_rotL_low,c2_rotL_low);

% FIGURE 5A --------------------------------------------------------------
% Input data
data1 = c1_rotE_low;
data2 = c2_rotE_low;
data3 = c3_rotE_low;
data = [data1;data2;data3];

limits = [min(data) max(data)];

% Group colours
cl1 = [222 177 97]./255;
cl2 = [53 155 124]./255;
cl3 = [34 102 141]./255;

% Plot
figure;

[a,b] = ksdensity(data1);

wdth = 0.5; % width of boxplot
% TODO, should probably be some percentage of max.height of kernel density plot

% density plot
d1 = area(b,a);
set(d1, 'FaceColor', cl1);
set(d1, 'EdgeColor', 'black');
set(d1, 'LineWidth', 1.5);
alpha(d1, 0.8);

hold on
[c,d] = ksdensity(data2);
d2 = area(d,c); 
set(d2, 'FaceColor', cl2);
set(d2, 'EdgeColor', 'black');
set(d2, 'LineWidth', 1.5);
alpha(d2, 0.6);

[e,f] = ksdensity(data3);
d3 = area(f,e);
set(d3, 'FaceColor', cl3);
set(d3, 'EdgeColor', 'black');
set(d3, 'LineWidth', 1.5);
alpha(d3, 0.4);

% make some space under the density plot for the boxplot
yl = get(gca,'YLim');
set(gca,'YLim',[-4 yl(2)]);

% jitter for raindrops
jit = (rand(size(data1)) - 0.5) * wdth/1.5;

% info for making boxplot
Y1 = quantile(data1,[0 0.25 0.75 0.5 0.02 1]);

% 'box' of 'boxplot'
b1 = rectangle('Position',[Y1(2) -0.5-(wdth*0.5) Y1(3)-Y1(2) wdth]);
set(b1,'EdgeColor','black')
set(b1,'LineWidth',1.5);
set(b1,'FaceColor',cl1+0.1);
% could also set 'FaceColor' here as Micah does, but I prefer without

% mean line
l1 = line([Y1(4) Y1(4)],[-0.5-(wdth*0.5) -0.5-(wdth*0.5)+wdth],'col','black','LineWidth',3);

% whiskers
l2 = line([Y1(3) Y1(6)], [-0.5 -0.5],'col','black','LineWidth',1.5);
l3 = line([Y1(1) Y1(2)],[-0.5 -0.5],'col','black','LineWidth',1.5);

% Group 2 boxplot
Y2 = quantile(data2,[0 0.25 0.75 0.5 0.02 1]);
b2 = rectangle('Position',[Y2(2) -1.15-(wdth*0.5) Y2(3)-Y2(2) wdth]);
set(b2,'EdgeColor','black')
set(b2,'LineWidth',1.5);
set(b2,'FaceColor',cl2+0.2);
% mean line
l4 = line([Y2(4) Y2(4)],[-1.15-(wdth*0.5) -1.15-(wdth*0.5)+wdth],'col','black','LineWidth',3);
% whiskers
l5 = line([Y2(3) Y2(6)], [-1.15 -1.15],'col','black','LineWidth',1.5);
l6 = line([Y2(1) Y2(2)],[-1.15 -1.15],'col','black','LineWidth',1.5);

% Group 3 boxplot
Y3 = quantile(data3,[0 0.25 0.75 0.5 0.02 1]);
b3 = rectangle('Position',[Y3(2) -1.8-(wdth*0.5) Y3(3)-Y3(2) wdth]);
set(b3,'EdgeColor','black')
set(b3,'LineWidth',1.5);
set(b3,'FaceColor',cl3+0.2);
% mean line
l7 = line([Y3(4) Y3(4)],[-1.8-(wdth*0.5) -1.8-(wdth*0.5)+wdth],'col','black','LineWidth',3);
% whiskers
l8 = line([Y3(3) Y3(6)], [-1.8 -1.8],'col','black','LineWidth',1.5);
l9 = line([Y3(1) Y3(2)],[-1.8 -1.8],'col','black','LineWidth',1.5);

% raindrops
s1 = scatter(data1,jit - 2.5);
s1.SizeData = 5;
s1.MarkerFaceColor = cl1+0.1;
s1.MarkerEdgeColor = 'none';
alpha(s1,0.05);

% Group 2
s2 = scatter(data2,jit - 3);
s2.SizeData = 5;
s2.MarkerFaceColor = cl2+0.2;
s2.MarkerEdgeColor = 'none';
alpha(s2,0.05);

% Group 3
s3 = scatter(data3,jit - 3.5);
s3.SizeData = 5;
s3.MarkerFaceColor = cl3+0.2;
s3.MarkerEdgeColor = 'none';
alpha(s3,0.05);

xlim([-(max(data)+0.1) max(data)+0.1]);
set(gca,'box','off','YTick',[],'YTickLabel',[],'FontSize',24,'FontName','Arial','linew',1.5);


% FIGURE 5B --------------------------------------------------------------
% Input data
data1 = c1_rotL_low;
data2 = c2_rotL_low;
data3 = c3_rotL_low;
data = [data1;data2;data3];

limits = [min(data) max(data)];

% Group colours
cl1 = [222 177 97]./255;
cl2 = [53 155 124]./255;
cl3 = [34 102 141]./255;

% Plot
figure;

[a,b] = ksdensity(data1);

wdth = 0.5; % width of boxplot
% TODO, should probably be some percentage of max.height of kernel density plot

% density plot
d1 = area(b,a);
set(d1, 'FaceColor', cl1);
set(d1, 'EdgeColor', 'black');
set(d1, 'LineWidth', 1.5);
alpha(d1, 0.8);

hold on
[c,d] = ksdensity(data2);
d2 = area(d,c); 
set(d2, 'FaceColor', cl2);
set(d2, 'EdgeColor', 'black');
set(d2, 'LineWidth', 1.5);
alpha(d2, 0.6);

[e,f] = ksdensity(data3);
d3 = area(f,e);
set(d3, 'FaceColor', cl3);
set(d3, 'EdgeColor', 'black');
set(d3, 'LineWidth', 1.5);
alpha(d3, 0.4);

% make some space under the density plot for the boxplot
yl = get(gca,'YLim');
set(gca,'YLim',[-4 yl(2)]);

% jitter for raindrops
jit = (rand(size(data1)) - 0.5) * wdth/1.5;

% info for making boxplot
Y1 = quantile(data1,[0 0.25 0.75 0.5 0.02 1]);

% 'box' of 'boxplot'
b1 = rectangle('Position',[Y1(2) -0.5-(wdth*0.5) Y1(3)-Y1(2) wdth]);
set(b1,'EdgeColor','black')
set(b1,'LineWidth',1.5);
set(b1,'FaceColor',cl1+0.1);
% could also set 'FaceColor' here as Micah does, but I prefer without

% mean line
l1 = line([Y1(4) Y1(4)],[-0.5-(wdth*0.5) -0.5-(wdth*0.5)+wdth],'col','black','LineWidth',3);

% whiskers
l2 = line([Y1(3) Y1(6)], [-0.5 -0.5],'col','black','LineWidth',1.5);
l3 = line([Y1(1) Y1(2)],[-0.5 -0.5],'col','black','LineWidth',1.5);

% Group 2 boxplot
Y2 = quantile(data2,[0 0.25 0.75 0.5 0.02 1]);
b2 = rectangle('Position',[Y2(2) -1.15-(wdth*0.5) Y2(3)-Y2(2) wdth]);
set(b2,'EdgeColor','black')
set(b2,'LineWidth',1.5);
set(b2,'FaceColor',cl2+0.2);
% mean line
l4 = line([Y2(4) Y2(4)],[-1.15-(wdth*0.5) -1.15-(wdth*0.5)+wdth],'col','black','LineWidth',3);
% whiskers
l5 = line([Y2(3) Y2(6)], [-1.15 -1.15],'col','black','LineWidth',1.5);
l6 = line([Y2(1) Y2(2)],[-1.15 -1.15],'col','black','LineWidth',1.5);

% Group 3 boxplot
Y3 = quantile(data3,[0 0.25 0.75 0.5 0.02 1]);
b3 = rectangle('Position',[Y3(2) -1.8-(wdth*0.5) Y3(3)-Y3(2) wdth]);
set(b3,'EdgeColor','black')
set(b3,'LineWidth',1.5);
set(b3,'FaceColor',cl3+0.2);
% mean line
l7 = line([Y3(4) Y3(4)],[-1.8-(wdth*0.5) -1.8-(wdth*0.5)+wdth],'col','black','LineWidth',3);
% whiskers
l8 = line([Y3(3) Y3(6)], [-1.8 -1.8],'col','black','LineWidth',1.5);
l9 = line([Y3(1) Y3(2)],[-1.8 -1.8],'col','black','LineWidth',1.5);

% raindrops
s1 = scatter(data1,jit - 2.5);
s1.SizeData = 5;
s1.MarkerFaceColor = cl1+0.1;
s1.MarkerEdgeColor = 'none';
alpha(s1,0.05);

% Group 2
s2 = scatter(data2,jit - 3);
s2.SizeData = 5;
s2.MarkerFaceColor = cl2+0.2;
s2.MarkerEdgeColor = 'none';
alpha(s2,0.05);

% Group 3
s3 = scatter(data3,jit - 3.5);
s3.SizeData = 5;
s3.MarkerFaceColor = cl3+0.2;
s3.MarkerEdgeColor = 'none';
alpha(s3,0.05);

xlim([-(max(data)+0.1) max(data)+0.1]);
set(gca,'box','off','YTick',[],'YTickLabel',[],'FontSize',24,'FontName','Arial','linew',1.5);


%% Energy Landscape

%create a bucket to store change in BOLD across all trials
MSDSig=[];

%loop across all trials and calc the change
%start and end of trials for analysis
% Change these parameters for different task conditions
strtTrial = 141;
endTrial = 220;

% Shuffle trial indices
trial_idx = [strtTrial:endTrial];
trial_perm = [];
for ii = 1:100
    trial_temp = trial_idx(randperm(size(trial_idx, 2)));
    trial_perm = [trial_perm; trial_temp];
end

% Define cluster group
% Change this parameter for different cluster
SUBJcluster = BIS_idx==3;

for ii =  1:size(trial_perm,1) % Iterate over all permutations

    MSD_temp = [];
    trial_idx2 = trial_perm(ii,:);

    for tt = 1:size(trial_idx2,2)
    
        for ss = 2:size(trial_idx2,2)

            disp([ii tt ss])
    
            trial1 = squeeze(sample_ts(trial_idx2(1),:,:));
            trial2 = squeeze(sample_ts(trial_idx2(ss),:,:));
    
            %displacement
            Disp = trial1-trial2;%order doesnt matter as square
            SqrdDisp = Disp.^2; %squaring
    
            %Optional if want to collape across regions
            MeanSqrdDisp = mean(SqrdDisp,2); % Average across regions (Time X 1)
    
            %just grab subjects of interest
            MeanSqrdDisp = MeanSqrdDisp(SUBJcluster);
    
    
            %store in bucket
            MSD_temp = [MSD_temp; MeanSqrdDisp(:)]; %dont care about regions/subjects concatenate all
    
        end
    
       trial_idx2 = circshift(trial_idx2,1,2);
    
    end

    MSDSig = [MSDSig MSD_temp];

end


%% Now we have msd vector we can calculate the energy (first pdf then NRG)

%clf
MSDSig = reshape(MSDSig,[],1);

% explore ds -- the vlaues of SD ---
figure;
histogram(MSDSig_E3,'Normalization','pdf') %how to choose
set(gca,'YScale','log');

%BM note test BW changing (equivalent to calculating to histogram)
% Adjust fit according to histogram plot
pd = fitdist(MSDSig_L3,'Kernel','BandWidth',0.1); % Fits probability distribution

% Change MSD range according to histogram plot
ds = 0.1:0.1:10;
pdfEstimate = pdf(pd,ds); % Returns probability density function

hold on
plot(ds,pdfEstimate)


%
nrg = -1.*log(pdfEstimate); % energy is negative log(pdf)

%

% plot(ds,nrg)
% xlabel('MSD change')
% ylabel('Energy')

%% Store nrg into separate variable for each task condition, cluster
nrg_L1 = nrg;
nrg_L2 = nrg;
nrg_L3 = nrg;


%% FIGURE 5C & D ---------------------------------------------------------

data1 = nrg_L1;
data2 = nrg_L2;
data3 = nrg_L3;

ds = 0.1:0.1:10;

% Group Colours
cl1 = [222 177 97]./255;
cl2 = [53 155 124]./255;
cl3 = [34 102 141]./255;

% Line plots
figure;
plot(ds,data1,'Color',cl1,'LineWidth',2.5);
hold on
plot(ds,data2,'Color',cl2,'LineWidth',2.5);
plot(ds,data3,'Color',cl3,'LineWidth',2.5);
% xlabel('MSD');
% ylabel('Energy')
set(gca,'box','off','FontSize',24,'FontName','Arial','LineWidth',1.5);
