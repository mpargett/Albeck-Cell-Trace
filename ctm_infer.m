%CTM_INFER
%   Univariate model inference for Cell Trace Models.  For data sets with
%   one independent and one dependent variable, ctm_infer calculates
%   information metrics between the variables and fits a set of models to
%   the data, with model discriminiation metrics by AIC (Akaike Information
%   Criterion).
%
%Usage:
%   W = ctm_infer(X, Y)
%       returns W, a structure containing all information and inference
%       results, as detailed in the 'Output' section below.   
%   W = ctm_infer(X, Y, 'PARAM1', VALUE, ...) uses the specified parameter
%       values (see 'Parameters' section below). 
%
%Inputs:
%   X - Independent variable data, as a column vector
%   Y - Dependent variable data, as a column vector
%
%Output:  
%   W.info  contains calculate information metrics, 'lcorr' for linear
%           correlation (Pearson), 'rcorr' for rank correlation (Spearman),
%           'mi' for mutual information, 'rmi' for relative mutual info
%           (scaled from 0 to 1).
%    .m     contains model fit results as a structure array with one
%           element per model type, and fields: 'name' the model name,
%           'model' the equation, 'sse' the fit sum of squared error, 
%           'np' the number of parameters, 'aicwt' the AIC weight against
%           the other models, 'evrat' the evidence ratio
%   .msum   summary function to print model fitness results concisely. Use
%           by typing 'W.msum(W.m)', replacing W with the name of your
%           output structure.
%
%Parameters:
%   RUNINFO - Logical flag to calculate information metrics, default = TRUE
%   RUNFITS - Logical flag to fit models, default = TRUE.
%   MODELS  - String or cell array of model names to fit, or name of a
%               designated library of models.  See ctm_modellib for model
%               and library names.  Default is 'general'.
%   

%Procedure:
%   Input handling
%   Calculate information metrics (univariate)
%   Fit linear models
%   Fit fixed non-linear models?
%   Compare models via AIC

function w = ctm_infer(x,y,varargin)
%Define default parameters
p.models = 'general';
p.mim = 'kNN_k';
p.runinfo = true;
p.runfits = true;

%Parse inputs
p = ct_input(varargin, p);

%Ensure column arrangement of inputs
sx = size(x); if length(sx) > 1 && sx(1) == 1; x = x(:); y = y(:); end

%Check for, and ensure, sorted data
if ~issorted(x); [x, xi] = sort(x, 'ascend'); y = y(xi); end

%% Information Metrics
if p.runinfo
%Linear correlation (Pearson)
[w.info.lcorr(1), w.info.lcorr(2)] = corr(x,y, 'type', 'Pearson', ...
    'rows', 'complete', 'tail', 'both');
%Rank correlation (Spearman)
[w.info.rcorr(1), w.info.rcorr(2)] = corr(x,y, 'type', 'Spearman', ...
    'rows', 'complete', 'tail', 'both');
%Mutual Information
%   Validate toolbox availability
if ~exist('HShannon_kNN_k_estimation', 'file'); 
    w.info.mi = [];     w.info.rmi = []; 
else
    %Calculate Shannon Entropies for each distribution and the joint
    %   Actually using Jaynes' Limiting Density of Discrete points
    %   adjustment: H(x) = - int( f(x)log(f(x)/m(x)) )dx, where m(x) is the
    %   sampling distribution on f(x).  Assuming a uniform sampling gives
    %   m(x) = 1/r, for r the range of x.  With Shannon's differential
    %   entropy as h(x), H(x) = h(x) - log(r).
    %       Note: Entropy from a uniform distribution is zero (0).
    %           All less entropic distributions have negative values.
    
    % 	IF want extra speed, use kNN version
    lrx = log(range(x)); lry = log(range(y));
    hest = eval(['HShannon_',p.mim,'_initialization(1);']); %#ok<NASGU>
    hxy  = eval(['HShannon_',p.mim,'_estimation([x,y]'', hest) - lrx - lry; ']);
    hx   = eval(['HShannon_',p.mim,'_estimation(x'', hest) - lrx;']);
    hy   = eval(['HShannon_',p.mim,'_estimation(y'', hest) - lry;']);
    %Caclulate mutual information (MI) and relative MI from entropies
    w.info.mi    = hx + hy - hxy;
    w.info.rmi   = w.info.mi./max(hx+lrx,hy+lry);
end
end


%% Build models
if p.runfits
%Robust univariate linear fit, by Theil-Sen
%   Use any separable univariate non-linear by Theil-Sen?
%   Exclude from AIC, of course?
[tsm, tsg] = ctm_theilsen(x,y,'w');
ts_w = struct('name', 'theil-sen', 'model', tsm, 'sse', tsg.sse, ...
    'np', 2);

%Get models to fit (using cti_modellib)
[mlib, mp] = ctm_modellib(p.models);
for s =1:size(mlib,1)  %FOR every model in the library
    %Get method to determine options to set
    fopt = fitoptions(mlib{s,2});   
    switch lower(fopt.Method)
        case 'linearleastsquares'
            [fo, gof] = fit(x, y, mlib{s,2});
        case 'nonlinearleastsquares'
            %Get start point estimates if available
            if ~isempty(mp{s}); mpar = mp{s}(x,y);
            else mpar = cell(1,3); end
            try
            [fo, gof] = fit(x, y, mlib{s,2}, 'StartPoint', mpar{3}, ...
                'Lower', mpar{1}, 'Upper', mpar{2});
            catch ME; fo = fittype(); gof = struct('sse',Inf);
                warning(ME.message);
            end
    end
    
    %   Assign fits and outputs for each model
    w.m(s).name     = mlib{s,1};
    w.m(s).model    = fo;
    w.m(s).sse      = gof.sse;
    w.m(s).np       = numcoeffs(fo);
end;    w.m = [ts_w, w.m];  %Prepend the Theil-Sen model


%% Evaluate evidence ratios by AIC
%Calculate AICc values for all models evaluated
ns = nnz(~isnan(y));  %Number of valid samples in dataset
AICc = (2.*[w.m.np] + ns.*log([w.m.sse]./ns) ...
    + 2*[w.m.np].*([w.m.np] + 1)./(ns - [w.m.np] - 1))';
dAIC = AICc - min(AICc);                %delta AICc
Awt = exp(-dAIC/2)./sum(exp(-dAIC/2));  %Akaike weights
[mxAwt] = max(Awt);                     %Find Max Akaike Wt
ER = Awt./mxAwt;                        %Evidence ratios vs. best?
%   Assign weights and ratios to output
Awt = num2cell(Awt); ER = num2cell(ER);
[w.m(:).aicwt] = deal(Awt{:});
[w.m(:).evrat] = deal(ER{:});

w.msum = @(x)[[{'#'},num2cell(1:numel(x))]',{'Model',x.name}', ...
    {'SSE',x.sse}', {'np',x.np}', {'AICwt',x.aicwt}', {'EvRat',x.evrat}'];

end

end

% %Testing
% dest = HShannon_kNN_k_initialization(1);
% n = ceil(10.^(2:0.25:4));
% figure;
% for s = 1:numel(n)
%     x = rand(n(s),1)*10;    y = 3 + 0.2.*x + rand(size(x))*0.5;
%     lrx = log(range(x)); lry = log(range(y));
%     hx(s) = HShannon_kNN_k_estimation(x', dest) - lrx;
%     hy(s) = HShannon_kNN_k_estimation(y', dest) - lry;
%     hxy(s) = HShannon_kNN_k_estimation([x,y]', dest) - lrx -lry;
%     z = rand(5,5000)*10;
%     h(s) = HShannon_kNN_k_estimation(z,dest) - log(prod(range(z,2)));
%     subplot(3,ceil(numel(n)/3),s);
%     mi = hx(s)+hy(s)-hxy(s);
%     rmi = mi ./ max(hx(s)+lrx,hy(s)+lry);
%     scatter(x,y); title(sprintf(['MI: ', num2str(mi), ...
%         ', RMI: ', num2str(rmi)]));
% end




