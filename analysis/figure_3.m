clear
close all
clc

%% Import data

A = importdata('data_experiment_2__by_participant.csv');
A = A.data;
B = importdata('data_experiment_2__with_features.csv');
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
id_cond      = B(:,2);
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
% s1_estimates(log10(s1_estimates)>3 & quest_id==3)=NaN;



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

%% Figure 3B: Effect of treatment on collective error

wfig(8,6)
rawebar(ze2(id_cond==3),ze2(id_cond==4))
alpha(.2)
box off
set(gca,'xtick',[1,2])
xticklabels({'Numbers','Fermi'})
xlabel 'Instruction'
ylabel 'Normalized error'
ylim([0, 2.2]);



print -dpdf figure3B.pdf

% First and second analyses: collective error 
dfer = id_cond == 4;
tbl = table(ze2,dfer,id_quest);
mdl = 'ze2 ~ dfer + (1|id_quest)';
M = fitlme(tbl,mdl);
computeCohen_d(ze2(id_cond==3),ze2(id_cond==4))

[h,p]=ranksum(ze2(id_cond==3),ze2(id_cond==4));

%% Figure 3C: rate of Fermi words per condition

wfig(8,6)
rawebar(WL_fermi(id_cond==3),WL_fermi(id_cond==4))
alpha(.2)
box off
set(gca,'xtick',[1,2])
xticklabels({'Numbers','Fermi'})
xlabel 'Instruction'
ylabel 'Rate of Fermi words'
ylim([0, .5]);

computeCohen_d(WL_fermi(id_cond==3),WL_fermi(id_cond==4))

print -dpdf figure3C.pdf

% Third analysis: Rate of Fermi words 
dfer = id_cond == 4;
tbl = table(WL_fermi,dfer,id_quest);
mdl = 'WL_fermi ~ dfer + (1|id_quest)';
M = fitlme(tbl,mdl);
[h,p]=ranksum(WL_fermi(id_cond==3),WL_fermi(id_cond==4));

%% Figure 3D: rate of numbers per condition

wfig(8,6)
rawebar(WL_nums(id_cond==3),WL_nums(id_cond==4))
alpha(.2)
box off
set(gca,'xtick',[1,2])
xticklabels({'Numbers','Fermi'})
xlabel 'Instruction'
ylabel 'Rate of s1 Numbers'
ylim([0, .5]);

computeCohen_d(WL_nums(id_cond==3),WL_nums(id_cond==4))

print -dpdf figure3D.pdf

% Fourth analysis: Rate of numbers 
dfer = id_cond == 4;
tbl = table(WL_nums,dfer,id_quest);
mdl = 'WL_nums ~ dfer + (1|id_quest)';
M = fitlme(tbl,mdl);
[h,p]=ranksum(WL_nums(id_cond==3),WL_nums(id_cond==4));

%% Figure 3E: WE similarity for each condition

wfig(8,6)
rawebar(WE_75(id_cond==3),WE_75(id_cond==4))
alpha(.2)
box off
set(gca,'xtick',[1,2])
xticklabels({'Numbers','Fermi'})
xlabel 'Instruction'
ylabel 'Word Embeddings Similarity'
ylim([.12 .22]);

computeCohen_d(WE_75(id_cond==3),WE_75(id_cond==4))

print -dpdf figure3E.pdf

% Fifth analysis: WE similarity
dfer = id_cond == 4;
tbl = table(WE_75,dfer,id_quest);
mdl = 'WE_75 ~ dfer + (1|id_quest)';
M = fitlme(tbl,mdl);
[h,p]=ranksum(WE_75(id_cond==3),WE_75(id_cond==4));


%% Effect of treatment on individual error

num_groups = unique(id_group(id_cond==3));
fer_groups = unique(id_group(id_cond==4));

cond_id = nan(size(discussed_id));
cond_id(ismember(group_id,num_groups))=3;
cond_id(ismember(group_id,fer_groups))=4;

num_dis = cond_id==3 & discussed_id;
num_und = cond_id==3 & ~discussed_id;
fer_dis = cond_id==4 & discussed_id;
fer_und = cond_id==4 & ~discussed_id;

% Error reduction on each condition

rnu = (ze1(num_und)-ze3(num_und));
rnd = (ze1(num_dis)-ze3(num_dis));
rfu = (ze1(fer_und)-ze3(fer_und));
rfd = (ze1(fer_dis)-ze3(fer_dis));


mrnu = nanmean(rnu);
mrnd = nanmean(rnd);
mrfu = nanmean(rfu);
mrfd = nanmean(rfd);

e_rnu = nanstd(rnu)/sqrt(sum(num_und));
e_rnd = nanstd(rnd)/sqrt(sum(num_dis));
e_rfu = nanstd(rfu)/sqrt(sum(fer_und));
e_rfd = nanstd(rfd)/sqrt(sum(fer_dis));



wfig(16,6)
subplot(1,2,1)
rawebar(rnu,rfu)
alpha(.2)
box off
set(gca,'xtick',[1,2])
xticklabels({'Numbers','Fermi'})
xlabel 'Instruction'
ylabel 'Error reduction'
ylim([-.25 2.5]);

subplot(1,2,2)
rawebar(rnd,rfd)
alpha(.2)
box off
set(gca,'xtick',[1,2])
xticklabels({'Numbers','Fermi'})
xlabel 'Instruction'
ylabel 'Error reduction'
ylim([-.25 2.5]);

print -dpdf figure3F.pdf



figure('color','w');
xticks = [1,1.25,3,3.25];
h1=plot(xticks(1)*[1,1],mrnu*[1,1],'ks','markerfacecolor','k','markersize',10);hold on;
plot(xticks(3)*[1,1],mrnd*[1,1],'ks','markerfacecolor','k','markersize',10);
h2=plot(xticks(2)*[1,1],mrfu*[1,1],'kd','markerfacecolor','k','markersize',10);
plot(xticks(4)*[1,1],mrfd*[1,1],'kd','markerfacecolor','k','markersize',10);
plot([0 5],[0,0],'k-')
plot(xticks(1)*[1,1],mrnu*[1,1]+e_rnu*[1,-1],'k-','linewidth',2);
plot(xticks(3)*[1,1],mrnd*[1,1]+e_rnd*[1,-1],'k-','linewidth',2);
plot(xticks(2)*[1,1],mrfu*[1,1]+e_rfu*[1,-1],'k-','linewidth',2);
plot(xticks(4)*[1,1],mrfd*[1,1]+e_rfd*[1,-1],'k-','linewidth',2);
ylim([-.1 .7])
xlim([xticks(1)-1,xticks(4)+1]);
set(gca,'xtick',[mean(xticks([1,2])),mean(xticks([3,4]))],'ytick',[0:.5:1.5],...
    'xticklabel',{'non-discussed','discussed'});
ylabel 'error reduction';
box off
legend([h1,h2],{'Numbers','Fermi'},'location','northwest')

%% STATS on Error Reduction

[h,p]=ranksum(rnd,rfd) 

er = ze1 - ze3;
dfer = cond_id == 4;
tbl = table(er,dfer,discussed_id,quest_id);
mdl = 'er ~ dfer*discussed_id + (1|quest_id)';
M = fitlme(tbl,mdl);


%% Use of Fermi Words (Whitelist) on each condition

mean(WL_nums(id_cond==3))
mean(WL_nums(id_cond==4))

[p,h]=ranksum(WL_fermi(id_cond==3),WL_fermi(id_cond==4))
[p,h]=ranksum(WL_nums(id_cond==3),WL_nums(id_cond==4))















