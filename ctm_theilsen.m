%CTM_THEILSEN
%   Theil-Sen linear model estimator, robust to outliers. 
%
%Usage:
%   [CF, GF] = ctm_theilsen(X, Y, METHOD, 'PARAM1', VALUE1, ...)
%
%Inputs:
%   X - Independent variable data, as a column vector
%   Y - Dependent variable data, as a column vector
%   METHOD - String indicating the averaging method to apply, with allowed
%       abbreviations in braces {}.
%       MEDIANALL{A}:      Takes the median across all pairwise slopes.
%       MEDIANMEDIAN{M}:   First takes the median of all slopes associated
%           with each data data point individually, then takes the median
%           of those medians.
%       WEIGHTEDMEDIAN{W}: Weights the use of each slope either using the
%           weight vector provided as a parameter, or (default) the
%           euclidean distance between points.
%
%Output:  
%   CF - fit model, as a cfit object.
%   GF - Goodness of fit metrics, including 
%       sse:        Sum of squared error
%       rsquare:    Coefficient of determination (R^2 or 'r squared')
%       dfe:        Statistical degrees of freedom (nSamples - nParameters)
%       adjrsquare: R^2 adjusted for degrees of freedom
%       rmse:       Root mean squared error
%
%Parameters:
%   WEIGHTS - Vector of weights to apply to residuals  during fitting. Will
%       preferentially weight fitness to some data. Must be size(X); the
%       weight applied to each slope is the product of the weights of the
%       two points involved.
%  

function [cf, gf] = ctm_theilsen(X,Y,method,varargin)
%Default parameters
p.weights = [];
%Parse additional inputs (parameter/value pairs)
p = ct_input(varargin, p);

%Assemble pairwise slopes
xys = pdist(Y(:),@minus)./pdist(X(:),@minus);

%Take median slope, per method requested
switch lower(method)
    case {'medianall', 'all', 'a'}
        m = median(xys,'omitnan');
    case {'medianmedian', 'mm', 'm'}
        m = median(median(squareform(xys),2,'omitnan'),1);
    case {'weightedmedian', 'wm', 'w'}
        if isempty(p.weights)   %IF weights are not provided
            [~,idx] = sort(X,'ascend');  X = X(idx);  Y = Y(idx);
            %Weight by inverse of local variance
            wt = pdist([X(:),Y(:)],'euclidean');
        else	wt = pdist(p.weights(:), @prod); %ELSE get pairwise weights
        end
        %Calculate weighted median
        [xys,idx] = sort(xys,'ascend');  wt = wt(idx);  %Sort slopes first
        ws = cumsum(wt);  m = xys(find(ws >= max(ws)/2, 1,'first'));
    otherwise; error('CTM_THEILSEN:BADMETHOD', ['Specified method is ',...
            'not valid.  Type "help ctm_theilsen" for options.']);
end

%Calculate intercept
b = median(Y - m*X,'omitnan');

%   Store X data for confidence estimate
nx = numel(X);  J = [X(:), ones(nx,1)];

%Calculate goodness of fit metrics
gf.sse = sum((Y - (m*X + b)).^2);
gf.rsquare = 1 - gf.sse/sum((Y - mean(Y,'omitnan')).^2);
gf.dfe = nx-2;
gf.adjrsquare = 1 - (1-gf.rsquare)*(nx-1)/(gf.dfe);
gf.rmse = sqrt( gf.sse/gf.dfe );

%Construct cfit object to encapsulate model
cf = cfit(fittype('poly1'), m, b, 'sse', gf.sse, 'dfe', gf.dfe, ...
    'Jacobian', J, 'xlim', [min(X), max(X)]);

end
