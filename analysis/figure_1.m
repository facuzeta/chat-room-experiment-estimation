clear
close all
clc

%% Import data

A = importdata('data_experiment_1__by_participant.csv');
A = A.data;
B = importdata('data_experiment_1__with_features.csv');
B = B.data;

%% Define variables

% data from individuals
participant_id  = A(:,1);
group_id        = A(:,2);
quest_id        = A(:,3);
s1_estimates    = A(:,4);
s1_conf         = A(:,5);
s3_estimates    = A(:,6);
s3_conf         = A(:,7);
n_num_s1        = A(:,8);
n_fer_wht       = A(:,9);
n_interv        = A(:,10);
n_words         = A(:,11);
discussed_id    = ~isnan(n_words);

% data from groups

id_group     = B(:,1);
id_quest     = B(:,3);
s1_mean      = B(:,5);
s2_estimates = B(:,6);
raters_fermi = B(:,7);
raters_nums  = B(:,8);
WL_fermi     = B(:,9);
WL_nums      = B(:,10);
WE_mean      = B(:,11);
WE_med       = B(:,15);
WE_75        = B(:,16);
WE_max       = B(:,17);
WE_90        = B(:,18);


% Correct answers
cor   = [171, 206, 373, 10800, 37, 7, 236, 48];

%% Normalize values

% In question #3, remove people who used a different scale (millions)
s1_estimates(log10(s1_estimates)>3 & quest_id==3)=NaN;



Nq = 8;
% predefine individual variables
z1   = nan(size(s1_estimates));
z3   = nan(size(s1_estimates));
ze1  = nan(size(s1_estimates));
ze3  = nan(size(s1_estimates));
e1   = nan(size(s1_estimates));
e3   = nan(size(s1_estimates));
% predefine group variables
z2   = nan(size(s2_estimates));
z1m  = nan(size(s2_estimates));
ze2  = nan(size(s2_estimates));
ze1m = nan(size(s2_estimates));
e2   = nan(size(s2_estimates));
e1m  = nan(size(s2_estimates));
% predefine correct numbers
cc = nan(Nq,1);
for q=1:Nq
    % normalize individual data
    indq = quest_id == q;
    medq = nanmedian(s1_estimates(indq));
    madq = mad(s1_estimates(indq),1);
    % normalize individual values and correct answers
    z1(indq)=(s1_estimates(indq)-medq)/madq;
    z3(indq)=(s3_estimates(indq)-medq)/madq;
    cc(q)   = (cor(q)-medq)./madq;
    % compute individual errors
    ze1(indq)=abs(z1(indq)-cc(q));
    e1(indq)=(abs(s1_estimates(indq)-cor(q)));
    ze3(indq)=abs(z3(indq)-cc(q));
    e3(indq)=(abs(s3_estimates(indq)-cor(q)));
    % normalize group data
    qind = id_quest == q;
    % normalize group values
    z2(qind)=(s2_estimates(qind)-medq)/madq;
    z1m(qind)=(s1_mean(qind)-medq)/madq;
    % compute group errors
    ze2(qind)=abs(z2(qind)-cc(q));
    e2(qind)=(abs(s2_estimates(qind)-cor(q)));
    ze1m(qind)=abs(z1m(qind)-cc(q));
    e1m(qind)=(abs(s1_mean(qind)-cor(q)));
    
end

% Remove outliers
thr = 2.5;
indout =  abs(z1)>thr | abs(z3)>thr;
z1(indout)=NaN;
ze1(indout)=NaN;
e1(indout)=NaN;
z3(indout)=NaN;
ze3(indout)=NaN;
e3(indout)=NaN;
indout2 =  abs(z2)>thr | abs(z1m)>thr;
z2(indout2)=NaN;
ze2(indout2)=NaN;
e2(indout2)=NaN;
z1m(indout2)=NaN;
ze1m(indout2)=NaN;
e1m(indout2)=NaN;


%% Figure 1A: Replication of Navajas et al., 2018
crowd_sizes = [4:4:32];
Nc = length(crowd_sizes);
Nit = 500;
ecrowd_z1 = nan(Nc,Nq,Nit);
ecrowd_z2 = nan(Nc,Nq,Nit);
ecrowd_z3 = nan(Nc,Nq,Nit);

for q=1:Nq
    z1q = z1(quest_id ==q & discussed_id & ~isnan(z1) & ~isnan(z3));
    z3q = z3(quest_id ==q & discussed_id & ~isnan(z1) & ~isnan(z3));
    z2q = z2(id_quest ==q & ~isnan(z2));
    for ic=1:Nc
        for it=1:Nit
            rand_ind = randperm(size(z1q,1));
            crowd_estimate_s1 = mean(z1q(rand_ind(1:crowd_sizes(ic))));
            crowd_estimate_s3 = mean(z3q(rand_ind(1:crowd_sizes(ic))));
            ecrowd_z1(ic,q,it) = abs(crowd_estimate_s1-cc(q));
            ecrowd_z3(ic,q,it) = abs(crowd_estimate_s3-cc(q));
            ind_rand = randperm(size(z2q,1));
            crowd_estimate_s2 = mean(z2q(ind_rand(1:crowd_sizes(ic)/4)));
            ecrowd_z2(ic,q,it) = abs(crowd_estimate_s2-cc(q));
        end
    end
end


m_ec_1 = mean(mean(ecrowd_z1,2),3);
m_ec_2 = mean(mean(ecrowd_z2,2),3);
m_ec_3 = mean(mean(ecrowd_z3,2),3);

e_ec_1 = std(mean(ecrowd_z1,2),[],3);
e_ec_2 = std(mean(ecrowd_z2,2),[],3);
e_ec_3 = std(mean(ecrowd_z3,2),[],3);

figure('color','w');
plot(crowd_sizes,m_ec_1,'bo-','markerfacecolor','b');
hold on;
plot(crowd_sizes,m_ec_2,'ko-','markerfacecolor','k');
plot(crowd_sizes,m_ec_1+e_ec_1,'b-');
plot(crowd_sizes,m_ec_1-e_ec_1,'b-');
plot(crowd_sizes,m_ec_2+e_ec_2,'k-');
plot(crowd_sizes,m_ec_2-e_ec_2,'k-');
xlim([0 35])
ylim([.25 1.35])
set(gca,'xtick',crowd_sizes,'ytick',[.4:.4:1.2])
xlabel 'crowd size'
ylabel 'normalized error'

%% Figure 1B: Effect of discussion on individual error
ze1d = nan(size(group_id,1),Nq);
ze3d = nan(size(group_id,1),Nq);
ze1u = nan(size(group_id,1),Nq);
ze3u = nan(size(group_id,1),Nq);
all_groups = unique(id_group);


for q=1:Nq
    d_groups = id_group(id_quest==q);
    u_groups = setdiff(all_groups,d_groups);
    these_indexes = ismember(group_id,d_groups) & quest_id==q;
    ze1d(these_indexes,q) = ze1(these_indexes);
    ze3d(these_indexes,q) = ze3(these_indexes);
    these_indexes = ismember(group_id,u_groups) & quest_id==q;
    ze1u(these_indexes,q) = ze1(these_indexes);
    ze3u(these_indexes,q) = ze3(these_indexes);
end

mze1d = nanmean(nanmean(ze1d,1));
mze3d = nanmean(nanmean(ze3d,1));
mze1u = nanmean(nanmean(ze1u,1));
mze3u = nanmean(nanmean(ze3u,1));

eze1d = nanstd(nanmean(ze1d,2));%/sqrt(size(ze1d,1));
eze3d = nanstd(nanmean(ze1d,2));%/sqrt(size(ze1d,1));
eze1u = nanstd(nanmean(ze1d,2));%/sqrt(size(ze1d,1));
eze3u = nanstd(nanmean(ze1d,2));%/sqrt(size(ze1d,1));
eze1d = mad(nanmean(ze1d,1),1);%/sqrt(size(ze1d,1));
eze3d = mad(nanmean(ze1d,1),1);%/sqrt(size(ze1d,1));
eze1u = mad(nanmean(ze1d,1),1);%/sqrt(size(ze1d,1));
eze3u = mad(nanmean(ze1d,1),1);%/sqrt(size(ze1d,1));

figure('color','w');
xticks = [1,1.25,3,3.25];
h1=plot(xticks([1,3]),[mze1d,mze3d],'k-','linewidth',2);
hold on;
h2=plot(xticks([2,4]),[mze1u,mze3u],'k--','linewidth',2);
plot(xticks([1,2]),[mze1d,mze1u],'bo','markerfacecolor','b');
plot(xticks([3,4]),[mze3d,mze3u],'ro','markerfacecolor','r');
plot(xticks(1)*[1,1],mze1d+[eze1d,-eze1d],'b-','linewidth',2);
plot(xticks(2)*[1,1],mze1u+[eze1u,-eze1u],'b-','linewidth',2);
plot(xticks(3)*[1,1],mze3d+[eze3d,-eze3d],'r-','linewidth',2);
plot(xticks(4)*[1,1],mze3u+[eze3u,-eze3u],'r-','linewidth',2);
xlim([xticks(1)-1,xticks(4)+1]);
ylim([0 1.75]);
legend([h1,h2],{'discussed','undiscussed'},'location','southwest')
set(gca,'xtick',[mean(xticks([1,2])),mean(xticks([3,4]))],'ytick',[0:.5:1.5],...
    'xticklabel',{'stage 1','stage 3'});
ylabel 'individual error';
box off

%% Figure 1C: collective estimates are better than the mean estimate

groups = unique(id_group);
Ngroups = length(groups);
for g=1:Ngroups
    mze2(g)=nanmean(ze2(id_group==groups(g)));
    mze1m(g)=nanmean(ze1m(id_group==groups(g)));
end

grey = .6*[1,1,1];
darkgrey = .2*[1,1,1];
figure('color','w');
plot(mze2,mze1m,'o','markersize',8,...
    'color',darkgrey,'markerfacecolor',grey);
hold on;
xlim([-0.25 2.75])
ylim([-0.25 2.75])
plot(xlim,xlim,'k--');
box off
set(gca,'xtick',[0:.5:2.5],'ytick',[0:.5:2.5]);
ylabel 'error of the group mean';
xlabel 'error of the group estimate';

figure('color','w')
[y,x]=hist(mze2-mze1m);
bar(x,y/sum(y),'facecolor',grey)
xlim([-2.5,2.5])
hold on
plot([0,0],ylim,'k--')
axis off

%% Figure 1D: influence of each group member as a function of initial accuracy

groups = unique(id_group);
Ngroups = length(groups);
z1c = nan(size(group_id,1),4);
z2c = nan(size(group_id,1),1);
j = 0;
for g=1:Ngroups
    for q=1:Nq
        these_indexes = group_id==groups(g) & quest_id==q;
        [~,i]=sort(ze1(these_indexes));
        aux=z1(these_indexes);
        j = j+1;
        if length(i)==4
            z1c(j,:) = aux(i);
        end
        my_index = id_group==groups(g) & id_quest==q;
        if sum(my_index)>0
            z2c(j)=z2(my_index);
        end
    end
end

ind_in=~isnan(sum([z1c,z2c],2));
z1c=z1c(ind_in,:);
z2c=z2c(ind_in);

M = fitlm(zscore(z1c,[],1),z2c);
betas = M.Coefficients.Estimate(2:5);
betas = betas./sum(betas);
b_SE = M.Coefficients.SE(2:5);

figure('color','w');
plot([4:-1:1],betas,'ko','markerfacecolor','k');
hold on
for k=1:4
plot(5-k*[1,1],betas(k)+b_SE(k)*[1,-1],'k-');
end
xlim([0 5])
% ylim([0 .7])
plot(xlim,.25*[1,1],'k--');
box off
set(gca,'xtick',[1:4],'ytick',[0:.2:.6],...
    'xticklabel',[])
xlabel 'individuals sorted by accuracy'
ylabel 'weight in group estimate'

