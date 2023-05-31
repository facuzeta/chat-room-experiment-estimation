function rawebar(data_y1,data_y2)

% Palette
c1 = [117,112,179]/255;
c2 = [27,158,119]/255;

% Calculate the mean and standard error of the mean for y1 and y2
mean_y1 = nanmean(data_y1);
sem_y1 = nanstd(data_y1) / sqrt(length(data_y1));
mean_y2 = nanmean(data_y2);
sem_y2 = nanstd(data_y2) / sqrt(length(data_y2));

% Define jitter amount
jitter_amount = 0.025;

% Add jitter to x-values
x1_jittered = ones(size(data_y1)) + jitter_amount*randn(size(data_y1));
x2_jittered = 2*ones(size(data_y2)) + jitter_amount*randn(size(data_y2));

s1=scatter(x1_jittered, data_y1,200,'Marker', '.', 'MarkerEdgeColor', c1); % Overlay y1 raw data with red dots and jittered x-values
hold on;
setMarkerColor(s1,c1,.2);
s2=scatter(x2_jittered, data_y2,200, 'Marker', '.', 'MarkerEdgeColor', c2); % Overlay y2 raw data with blue dots and jittered x-values
setMarkerColor(s2,c2,.2);
plot([1,1],mean_y1+sem_y1*[1,-1],'-','color','k','LineWidth',1);
% plot([1,1],mean_y1,'o','color','k','markerfacecolor',c1,'markersize',10);
plot([1,1]+.05*[1,-1],mean_y1*[1,1],'-','color','k','LineWidth',2);
plot([2,2],mean_y2+sem_y2*[1,-1],'-','color','k','LineWidth',1);
% plot([2,2],mean_y2,'o','color','k','markerfacecolor',c2,'markersize',10);
plot([2,2]+.05*[1,-1],mean_y2*[1,1],'-','color','k','LineWidth',2);

xlim([0.5, 2.5]);

end