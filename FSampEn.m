function result = FSampEn(data, r, dim)

%% Fast Sample Entropy calculation

% data = EEG time series data [EEG.data]
% r = tolerance => percent of std to accept as a match
% dim = embedding dimension => length of comparison window

% create matrix of data by embedding dimension
data_matrix = NaN(dim+1, length(data));
for a = 1:dim+1
    data_matrix(a,1:length(data)+1-a) = data(a:end);
end
data_matrix = data_matrix';

% take pairwise distances across each window
matrix_B = pdist(data_matrix(:,1:dim), 'chebychev');
matrix_A = pdist(data_matrix(:,1:dim+1), 'chebychev');

% take sum of counts of distances that fall within tolerance range
B = sum(matrix_B <= r*std(data));
A = sum(matrix_A <= r*std(data));

% calculate ratio of natural logarithm, with normalization correction
ratio = -log((A/B))*((length(data)-dim+1)/(length(data)-dim-1));

% corrections
max = -log(2/((length(data)-dim-1)*(length(data)-dim)));
if isinf(ratio)
    ratio = max;
end
if(ratio>max)
    ratio = max;
end

result = ratio;
end