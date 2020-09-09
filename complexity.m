%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes various measures of complexity for a %
% time series.                                                %
%                                                             %
% INPUTS                                                      %
%  Parameters:                                                %
%   series = time series                                      %
%   type = string denoting type of complexity measure:        %
%        'SE' = Sample Entropy                                %
%        'HHSE' = Hilbert-Huang Spectral Entropy              %
%        'HFD' = Higuchi Fractal Dimension                    %
%  Optional:                                                  %
%   dim = dimensional modifier for SE and HFD                 %
%   r = tolerance modifier for SE                             %
%   fs = sample rate for HHSE                                 %
%                                                             %
% OUTPUTS                                                     %
%   result = single complexity measure value                  %
%                                                             %
%                                                             %
%                 Created by Brian Lord                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = complexity(series,type, varargin)
%% Parameters
% defaults:
r = 0.2;
dim = 3;
fs = 256;
data = series;
%data = (data-mean(data))/std(data); 

if nargin < 2
    error('Missing parameters!');
elseif nargin == 2
    disp('Using default parameters...');
    if strcmp(type,'HFD')
        dim = 5;
    end
elseif nargin == 3
    dim = varargin{1};
elseif nargin == 4
    dim = varargin{1};
    r = varargin{2};
elseif nargin ==5
    dim = varargin{1};
    r = varargin{2};
    fs = varargin{3};
end
    
        

%% Sample Entropy %%
% The Sample Entropy is the negative average natural logarithm of the
% conditional probability that two sequences that are similar for m points
% remain similar within a tolerance r at the next point, with self-matches
% not included in the probability calculation. It is a time measure of
% complexity.
%
% Richman, J. S., & Moorman, J. R. (2000). Physiological time-series 
%         analysis using approximate entropy and sample entropy. Am J 
%         Physiol Heart Circ Physiol, 278, H2039-H2049. 
%         doi:10.1152/ajpheart.2000.278.6.H2039
% Keller, K., Unakafov, A., & Unakafova, V. (2014). Ordinal Patterns, 
%         Entropy, and EEG. Entropy, 16(12), 6212-6239. 
%         doi:10.3390/e16126212
%
% r = tolerance => percent of std to accept as a match (default = 0.2)
% dim = embedded dimension => length of comparison window (default = 3)

if strcmp(type,'SE')
    
    data = (data-mean(data))/std(data); % Normalize time series
    m = dim; % embedded dimension
    N = length(data);
    r = r*std(data);
    
    sumcounts = zeros(1,2);
    for a = m:m+1
        counts = zeros(1,(N-a)); % template for counting matches
        for j = 1:(N-m)
            tempWindow = data(j:(j+a-1)); % Data window to compare
            distances = zeros(1,(N-a)); % template for distances
            for k = 1:(N-m)
                tempWindow2 = data(k:(k+a-1)); % Sliding data window
                distances(k) = max(abs(tempWindow - tempWindow2)); % Compute chebychev distance
            end
            distances(j) = []; % Eliminate self-matches
            counts(j) = sum(distances<=r); % Count matches that were within tolerance
        end
        sumcounts(a-m+1) = sum(counts)*(2/((N-a-1)*(N-a))); % Compute
    end
    
    SampEn = log(sumcounts(1))-log(sumcounts(2));
    result = SampEn;
    
end


%% Hilbert-Huang Spectral Entropy
% The Hilbert-Huang Spectral Entropy applies the Shannon entropy to the
% Hilbert-Huang spectrum of the signal. It is a time-frequency measure of
% complexity.

% Huang, N. E. et al. (1998). The empirical mode decomposition and the 
%         Hilbert spectrum for nonlinear and non-stationary time series 
%         analysis. Proc R Soc Lond A, 454, 903-995. 
%         doi:10.1098/rspa.1998.0193
% Li, X. et al. (2008). Analysis of depth of anesthesia with Hilbert-Huang 
%         spectral entropy. Clin Neurophysiol, 119(11), 2465-2475. 
%         doi:10.1016/j.clinph.2008.08.006
%
% fs = sample rate (default = 256)

if strcmp(type,'HHSE')
    
    imf = emd(data,'Display',0); % Perform empirical mode decomposition
    % imf = intrinsic mode functions
    
    HS = full(hht(imf,fs)); % Get Hilbert spectrum
    nHS = HS/sum(HS); % Normalize HS
    nHS(nHS==0) = []; % remove zero numbers
    HHSE = -sum(nHS.*log10(nHS)); % Calculate Shannon entropy of spectrum
    result = HHSE;
    
end


%% Higuchi Fractal Dimension
% The Higuchi Fractal Dimension describes the degree of self-similarity
% along a 1D time series, according to a sequencing integer, k. The
% function returns the dimensionality value, D, where 1 <= D <= 2. A higher
% D corresponds to higher complexity. It is a time measure of complexity.
%
% Higuchi, T. (1988). Approach to an Irregular Time Series on the Basis 
%         of the Fractal Theory. Physica D, 31, 277-283.
% Accardo, A., Affinito, M., Carrozzi, M., & Bouquet, F. (1997). Use of 
%         the fractal dimension for the analysis of electroencephalographic
%         time series. Biol Cybern, 77, 339-350. 
%         doi:10.1007/s004220050394
%
% dim = kmax => maximum size of iteration (default = 5)

if strcmp(type,'HFD')
    
    % Initialize components
    kmax = dim;
    N = length(data);
    L = zeros(1,kmax);
    x = zeros(1,kmax);
    y = zeros(1,kmax);
    
    for k = 1:kmax
        for m = 1:k
            norm_factor = (N-1)/(round((N-m)/k)*k); % normalization factor
            X = sum(abs(diff(data(m:k:N)))); % generate iterative sum difference
            L(m)=X*norm_factor/k; % mean of normalized differences
        end

        y(k)=log(sum(L)/k); % ln(L(k))
        x(k)=log(1/k); % ln(1/k)
    end
    
    D = polyfit(x,y,1); % linear regression fit
    HFD = D(1); % return slope
    result = HFD;
   
end

end
