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



%% Figure 2B: Rating of implementing Fermi method and group error

figure('color','w','paperunits','centimeters','paperposition',[1 1 8 6]);
rs = prctile(raters_fermi,[0:20:100]);
m_ze2 = nan(length(rs)-1,1);
e_ze2 = nan(length(rs)-1,1);
nrs = nan(length(rs)-1,1);
xrs = nan(length(rs)-1,1);
e_xrs = nan(length(rs)-1,1);
for r=1:length(rs)-1
    indexes = raters_fermi>=rs(r) & raters_fermi<rs(r+1);
    m_ze2(r)=nanmean(ze2(indexes));
    e_ze2(r)=nanstd(ze2(indexes))/sqrt(sum(indexes));
    nrs(r)=sum(indexes);
    xrs(r)=nanmean(raters_fermi(indexes));
    e_xrs(r)=nanstd(raters_fermi(indexes))/sqrt(sum(indexes));
end
plot(xrs,m_ze2,'ko','markerfacecolor','k','markersize',10);hold on;
for ir=1:length(rs)-1
plot([xrs(ir)-e_xrs(ir),xrs(ir)+e_xrs(ir)],m_ze2(ir)*[1,1],'k-');
plot(xrs(ir)*[1,1],m_ze2(ir)*[1,1]+e_ze2(ir)*[-1,1],'k-');
end
mdl = fitlm(raters_fermi,ze2);
mdl2 = fitlm(zscore(raters_fermi),ze2);
xi = [0:.1:10];
[yi,yci]=predict(mdl,xi');
plot(xi,yi,'r-','linewidth',2);
plot(xi,yci(:,1),'r-');
plot(xi,yci(:,2),'r-');
ylim([0.35 1.25])
xlim([-0.25 10.25])
xlabel 'Mean rating - Fermi'
ylabel 'Normalized error'
set(gca,'xtick',[0:2:10],'ytick',[.4:.2:1.2],'xticklabel',[0:2:10])
box off
print -dpdf figure2B.pdf

%% Figure 2C: Rating of implementing Numbers method and group error

figure('color','w','paperunits','centimeters','paperposition',[1 1 8 6]);
rs = prctile(raters_nums,[0:20:100]);
m_ze2 = nan(length(rs)-1,1);
e_ze2 = nan(length(rs)-1,1);
nrs = nan(length(rs)-1,1);
xrs = nan(length(rs)-1,1);
e_xrs = nan(length(rs)-1,1);
for r=1:length(rs)-1
    indexes = raters_nums>=rs(r) & raters_nums<rs(r+1);
    m_ze2(r)=nanmean(ze2(indexes));
    e_ze2(r)=nanstd(ze2(indexes))/sqrt(sum(indexes));
    nrs(r)=sum(indexes);
    xrs(r)=nanmean(raters_nums(indexes));
    e_xrs(r)=nanstd(raters_nums(indexes))/sqrt(sum(indexes));
end
plot(xrs,m_ze2,'ko','markerfacecolor','k','markersize',10);hold on;
for ir=1:length(rs)-1
plot([xrs(ir)-e_xrs(ir),xrs(ir)+e_xrs(ir)],m_ze2(ir)*[1,1],'k-');
plot(xrs(ir)*[1,1],m_ze2(ir)*[1,1]+e_ze2(ir)*[-1,1],'k-');
end
mdl = fitlm(raters_nums,ze2);
mdl2 = fitlm(zscore(raters_nums),ze2);
xi = [0:.1:10];
[yi,yci]=predict(mdl,xi');
plot(xi,yi,'r-','linewidth',2);
plot(xi,yci(:,1),'r-');
plot(xi,yci(:,2),'r-');
ylim([0.35 1.25])
xlim([-0.25 10.25])
xlabel 'Mean rating - Numbers'
ylabel 'Normalized error'
set(gca,'xtick',[0:2:10],'ytick',[.4:.2:1.2],'xticklabel',[0:2:10])
box off
print -dpdf figure2C.pdf


%% Figure 2E: Whitelist of Fermi keywords and group error
ep = 10e-3;
zWL_fermi = zscore(norminv(WL_fermi+ep));
figure('color','w','paperunits','centimeters','paperposition',[1 1 8 6]);
rs = prctile(zWL_fermi,[0:20:100]);
m_ze2 = nan(length(rs)-1,1);
e_ze2 = nan(length(rs)-1,1);
nrs = nan(length(rs)-1,1);
xrs = nan(length(rs)-1,1);
e_xrs = nan(length(rs)-1,1);
for r=1:length(rs)-1
    indexes = zWL_fermi>=rs(r) & zWL_fermi<rs(r+1);
    m_ze2(r)=nanmean(ze2(indexes));
    e_ze2(r)=nanstd(ze2(indexes))/sqrt(sum(indexes));
    nrs(r)=sum(indexes);
    xrs(r)=nanmean(zWL_fermi(indexes));
    e_xrs(r)=nanstd(zWL_fermi(indexes))/sqrt(sum(indexes));
end
plot(xrs,m_ze2,'ko','markerfacecolor','k','markersize',10);hold on;
for ir=1:length(rs)-1
plot([xrs(ir)-e_xrs(ir),xrs(ir)+e_xrs(ir)],m_ze2(ir)*[1,1],'k-');
plot(xrs(ir)*[1,1],m_ze2(ir)*[1,1]+e_ze2(ir)*[-1,1],'k-');
end
mdl = fitlm(zWL_fermi,ze2);
mdl2 = fitlm(zscore(zWL_fermi),ze2);
xi = [-5:.1:5];
[yi,yci]=predict(mdl,xi');
plot(xi,yi,'r-','linewidth',2);
plot(xi,yci(:,1),'r-');
plot(xi,yci(:,2),'r-');
xlim([-2.5 2.5])
ylim([0.35 1.25])
xlabel 'Whitelist Fermi'
ylabel 'Normalized error'
set(gca,'xtick',[0:2:10],'ytick',[.4:.2:1.2])
box off

print -dpdf figure2E.pdf

%% Figure 2E: Whitelist of numbers and group error
ep = 10e-3;
zWL_nums = zscore(norminv(WL_nums+ep));
figure('color','w','paperunits','centimeters','paperposition',[1 1 8 6]);
rs = prctile(zWL_nums,[0:20:100]);
m_ze2 = nan(length(rs)-1,1);
e_ze2 = nan(length(rs)-1,1);
nrs = nan(length(rs)-1,1);
xrs = nan(length(rs)-1,1);
e_xrs = nan(length(rs)-1,1);
for r=1:length(rs)-1
    indexes = zWL_nums>=rs(r) & zWL_nums<rs(r+1);
    m_ze2(r)=nanmean(ze2(indexes));
    e_ze2(r)=nanstd(ze2(indexes))/sqrt(sum(indexes));
    nrs(r)=sum(indexes);
    xrs(r)=nanmean(zWL_nums(indexes));
    e_xrs(r)=nanstd(zWL_nums(indexes))/sqrt(sum(indexes));
end
plot(xrs,m_ze2,'ko','markerfacecolor','k','markersize',10);hold on;
for ir=1:length(rs)-1
plot([xrs(ir)-e_xrs(ir),xrs(ir)+e_xrs(ir)],m_ze2(ir)*[1,1],'k-');
plot(xrs(ir)*[1,1],m_ze2(ir)*[1,1]+e_ze2(ir)*[-1,1],'k-');
end
mdl = fitlm(zWL_nums,ze2);
mdl2 = fitlm(zscore(zWL_nums),ze2);
xi = [-5:.1:5];
[yi,yci]=predict(mdl,xi');
plot(xi,yi,'r-','linewidth',2);
plot(xi,yci(:,1),'r-');
plot(xi,yci(:,2),'r-');
xlim([-2.5 2.5])
ylim([0.35 1.25])
xlabel 'Whitelist Numbers'
ylabel 'Normalized error'
set(gca,'xtick',[0:2:10],'ytick',[.4:.2:1.2])
box off

print -dpdf figure2F.pdf

%% Figure 2H: Word embeddings similarity and group error
figure('color','w','paperunits','centimeters','paperposition',[1 1 8 6]);
rs = prctile(WE_75,[0:20:100]);
m_ze2 = nan(length(rs)-1,1);
e_ze2 = nan(length(rs)-1,1);
nrs = nan(length(rs)-1,1);
xrs = nan(length(rs)-1,1);
e_xrs = nan(length(rs)-1,1);
for r=1:length(rs)-1
    indexes = WE_75>=rs(r) & WE_75<rs(r+1);
    m_ze2(r)=nanmean(ze2(indexes));
    e_ze2(r)=nanstd(ze2(indexes))/sqrt(sum(indexes));
    nrs(r)=sum(indexes);
    xrs(r)=nanmean(WE_75(indexes));
    e_xrs(r)=nanstd(WE_75(indexes))/sqrt(sum(indexes));
end
plot(xrs,m_ze2,'ko','markerfacecolor','k','markersize',10);hold on;
for ir=1:length(rs)-1
plot([xrs(ir)-e_xrs(ir),xrs(ir)+e_xrs(ir)],m_ze2(ir)*[1,1],'k-');
plot(xrs(ir)*[1,1],m_ze2(ir)*[1,1]+e_ze2(ir)*[-1,1],'k-');
end
mdl = fitlm(WE_75,ze2);
mdl2 = fitlm(nanzscore(WE_75),ze2);
xi = [.13:.01:.21];
[yi,yci]=predict(mdl,xi');
plot(xi,yi,'r-','linewidth',2);
plot(xi,yci(:,1),'r-');
plot(xi,yci(:,2),'r-');
xlim([.13 .21])
ylim([0.35 1.25])
xlabel 'Word Embeddings Similarity'
ylabel 'Normalized error'
set(gca,'xtick',[0.2:.02:.35],'ytick',[.4:.2:1.2])
box off

print -dpdf figure2H.pdf


%% Figure 2I: Correlation matrix for different measures

Xdata = [raters_fermi,WL_fermi,WE_75,raters_nums,WL_nums];
[r,p]=corr(Xdata,'type','spearman');
labels = {{'rating';'Fermi'},{'whitelist';'Fermi'},{'word';'similarity'},{'rating';'Numbers'},{'whitelist';'Numbers'}};
r(r>.99)=NaN;

figure('color','w','paperunits','centimeters','paperposition',2*[1 1 8 6]);

% scatter plot
n = size(r, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
scatter(x(:), y(:), 2500.*abs(r(:)), r(:), 'filled', 'MarkerFaceAlpha', 0.6)


text(2:n, (n:-1:2) + 0.5, labels(1:end-1), 'HorizontalAlignment', 'right','FontAngle','italic')
text((2:n) + 0.5, repmat(n + 1, n-1, 1)+.3, labels(2:end), ...
    'HorizontalAlignment', 'center', 'Rotation', 0,'FontAngle','italic')
for k=1:n*n
if isnan(r(k))
str{k}='';
else
if p(k)<.05
str{k}=[num2str(round(r(k)*100)/100), '*'];
else
str{k}=[num2str(round(r(k)*100)/100)];
end
end
end
hold on
textscatter(x(:), y(:), str,'FontWeight','bold')


h = gca;
j=colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.9;
axis(h, 'equal')
j.Position = j.Position+[0 0.05 0 0.05];
j.Ticks = [-.4:.2:.4];
colormap(redblue)
j.Visible = 'off';
clim([-.5 .5])

print -dpdf figure2I.pdf

