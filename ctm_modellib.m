%CTI_MODELLIB
%   Library of non-linear models for consideration with Cell Trace data.
%
%Usage:
%   [M, MP] = ctm_modellib(TYPE)
%       returns W, a structure containing all information and inference
%       results, as detailed in the 'Output' section below.   
%
%Inputs:
%   TYPE - String or cell array list of model types or model libraries.
%
%Outputs:  
%   M   - Models, as a nModels x 2 cell array. Column 1 contains model
%           names, column 2 contains fittype objects.
%   MP  - Model fitting parameter functions, as a nModels x 1 cell array.
%           When called as in P = MP(X,Y) with X and Y being data for
%           fitting, each returns a 1 x 3 cell array of 
%           {LowerBounds, UpperBounds, StartPoint}.  Only used for
%           non-linear models.
%
%Available Models:
%   Any model type defined for MATLAB's FITTYPE.
%   EXPB, EXPB2 - One- and two-term exponential models with baseline.
%   EXPD, EXPD2 - One- and two-term exponential models with baseline and a
%                   delayed onset.
%   HILL, HILL2 - One- and two-term hill functions.
%
%Available Libraries:
%   GENERAL(default) -  Poly1-2, Fourier1-3, Gauss1-3, EXPB1, EXPD1, HILL1
%   SIMPLE      -  	Poly1-2, EXPB1, EXPD1, HILL1
%   EXPONENTIAL -   EXP1-2, EXPB1-2, EXPD1-2
%   SATURATING  -   EXPB1-2, EXPD1-2, HILL1-2
%   BASIS       -   Poly1-9, Fourier1-8, Gauss1-8
%   EXHAUSTIVE  -   Poly1-9, Fourier1-8, Gauss1-8, Power1-2, Exp1-2,
%                   EXPB1-2, EXPD1-2, HILL1-2
%
%NOTE: Fitting of EXPB2, EXPD2, and HILL2 untested

function [m, mp] = ctm_modellib(t)
%Assembles models requested by type input, t
if ~iscell(t); t = {t}; end

m = [];  mp = [];
for s = 1:numel(t)
    %Parse model type (is list, or is indiv)
    %   Check for matching list
    lst = list_def(t{s});
    %   If no list match, check model names and retreive def
    if isempty(lst); [mdl, mfun] = model_def(t{s}); 
    else [mdl, mfun] = model_def(lst);    end  %If matching list, retrieve models
    %Append to output
    m = [m; mdl];   mp = [mp; mfun];    %#ok<AGROW>
    
end

end


% ------------------------------------------------------------------------
%             ----- --- -----  Subfunctions  ----- --- -----
% ------------------------------------------------------------------------

%% List definitions
function lst = list_def(n)
lst = [];
switch lower(n)
    case 'general'
        lst = {'poly1', 'poly2', 'fourier1', 'fourier2', 'fourier3', ...
            'gauss1', 'gauss2', 'gauss3', 'expb1', 'expd1', 'hill1'};
    case 'simple'
        lst = {'poly1', 'poly2', 'expb1', 'expd1', 'hill1'};
    case 'exponential'
        lst = {'exp1', 'exp2', 'expb1', 'expb2', 'expd1', 'expd2'};
    case 'saturating'
        lst = {'expb1', 'expb2', 'expd1', 'expd2', 'hill1', 'hill2'};
    case 'basis'
        lst = [cellfun(@(x)['poly',num2str(x)],num2cell(1:9),'un',0),...
               cellfun(@(x)['fourier',num2str(x)],num2cell(1:8),'un',0),...
               cellfun(@(x)['gauss',num2str(x)],num2cell(1:8),'un',0)];
    case 'exhaustive'
        lst = [cellfun(@(x)['poly',num2str(x)],num2cell(1:9),'un',0),...
               cellfun(@(x)['fourier',num2str(x)],num2cell(1:8),'un',0),...
               cellfun(@(x)['gauss',num2str(x)],num2cell(1:8),'un',0),...
               {'power1', 'power2', 'exp1', 'exp2', 'expb1', 'expb2', ...
               'expd1', 'expd2', 'hill1', 'hill2'}];
end

end


%% Model definitions (names/handles for cfit, and customs)
function [mn, mp] = model_def(n)
if ~iscell(n); n = {n}; end     %Ensure cell encapsulation

%Match names
mn = cell(numel(n),2);  mp = cell(numel(n),1);
fo = fitoptions('Method','Nonlinearleastsquares', 'Display', 'off');
for s = 1:numel(n);    n{s} = lower(n{s});  mn{s,1} = n{s};
    %   Check if model name is a proper cfit name
    if ~isempty(regexpi(n{s}, ['poly[1-9]|fourier[1-8]|gauss[1-8]',...
            '|exp[1-2]|power[1-2]|rat\d{2}|sin[1-8]|weibull']))
        mn{s,2} = fittype(n{s}); 	%Set fittype object
    else
        switch lower(n{s})
            %Specify custom names and definitions
            case {'expb1'};       %Exponential (allowing baseline)
                mn{s,2} = fittype(@(c0,a,b,x)c0 + a*exp(b*x), ...
                    'coefficients', {'c0','a','b'}, 'options', fo);
            case {'expb2'};     %Exponential w/ 2 terms
                mn{s,2} = fittype(@(c0,a1,b1,a2,b2,x)c0 + ...
                    a1*exp(b1*x) + a2*exp(b2*x), 'coefficients', ...
                    {'c0','a1','b1','a2','b2'}, 'options', fo);
            case {'expd1'};      %Exponential with delay
                mn{s,2} = fittype(@(c0,a,b,t0,x)c0 + ...
                    a*(double(x<t0) + (x>=t0).*exp(b*(x-t0))), ...
                    'coefficients', {'c0','a','b','t0'}, 'options', fo);
            case {'expd2'};    %Exponential with delay, 2 terms
                mn{s,2} = fittype(@(c0,a1,b1,t01,a2,b2,t02,x)c0 + ...
                    a1*(double(x<t01) + (x>=t01).*exp(b1*(x-t01))) + ...
                    a2*(double(x<t02) + (x>=t02).*exp(b2*(x-t02))), ...
                    'coefficients', {'c0','a1','b1','t01','a2','b2','t02'},...
                    'options', fo);
            case {'hill1'};      %Hill equation
                mn{s,2} = fittype(@(c0,a,ka,h,x)c0 + a./(1 + (ka./x).^h), ...
                    'coefficients', {'c0','a','ka','h'}, 'options', fo);
            case {'hill2'};    %Hill equation w/ 2 terms
                mn{s,2} = fittype(@(c0,a1,ka1,h1,a2,ka2,h2,x)c0 + ...
                    a1./(1 + (ka1./x).^h1) + a2./(1 + (ka2./x).^h2), ...
                    'coefficients', {'c0','a1','ka1','h1','a2','ka2','h2'},...
                    'options', fo);
            otherwise;  warning('Bad Model Name');  continue;
        end
    end
    fopt = fitoptions(mn{s,2});
    if strcmpi(fopt.Method, 'nonlinearleastsquares')
        mp{s} = @(x,y)mdlstart(n{s},x,y);  %Assign start point fn
    end
end

end


%% Model start point estimators
function mp = mdlstart(typ, x, y)
%Parse model name
mn = regexpi(typ, '(?<m>^\D*)(?<n>\d?)', 'names');
%   Correct for model types w/o number, and get double values
if isempty(mn.n); mn.n = 1; else mn.n = str2double(mn.n); end    

%Range of output, with a safety factor
ry = range(y)*10;

%Use smoothed data for estimates
sy = smooth(x, y, 0.2, 'loess'); 
%Get derivatives and estimate if an exponential model
if any(strcmpi(mn.m, {'expb','expd'}))
    dx = diff(x);  dy = diff(sy)./dx;  dxb = 1./[dx(1:end-1),dx(2:end)];
    dy = sum([dy(1:end-1),dy(2:end)] .* ...
        bsxfun(@rdivide,dxb, sum(dxb,2)),2);
    %Identify sign of a and b terms
    ns = numel(sy);  nd = numel(dy);
    %   Sign of second derivative defines sign of a
    as = sign(median( dy(ceil(nd/2)+1:end) - dy(1:floor(nd/2)) ));
    %   Sign of first derivative is sign(a)*sign(b)
    bs = sign(median( sy(ceil(ns/2)+1:end) - sy(1:floor(ns/2)) )) * as;
    %Iteratively get 2-point estimate
    b = bs; be = 0;
    while abs(b-be) > 1e-3;     be = b; %Retain last new b
        a = (sy(end)-sy(1))./(exp(b*x(end))-exp(b*x(1)));       %Est. a
        if be == 0; c = sy(1) - a;                              %Init. c
        else        c = sy(1) - a*exp(b*x(1)); end;             %Est. c
        aa = log((sy-c)./a)./x;  b = median( aa(~imag(aa)) );   %Est. b
    end
    if any([sign([a,b])~=[as,bs], isnan(c)]);  %Catch for poor convergence
        b = bs; a = (sy(end)-sy(1))./(exp(b*x(end))-exp(b*x(1)));
        c = sy(1) - a;      %Uses first pass guesses
    end
end

bl = []; bu = []; sp = [];
switch mn.m
    case 'expb'      %Basic exponential (w/ baseline)
        %   Coefs: c0, a1, b1
        sp = [c,a,b];  %Set start
        sp = [sp(1), repmat(sp(2:end), 1, mn.n)];
        %Set coefficient bounds
        bl = [-ry, -ry, -Inf]; bl = [bl(1), repmat(bl(2:end), 1, mn.n)];
        bu = [ry,  ry,  Inf]; bu = [bu(1), repmat(bu(2:end), 1, mn.n)];
        %Set both to identical for Two-term exponential ??
    case 'expd'     %Exponential with delay
        %   Coefs: c0, a1, b1, td
        y1 = sy(1:ceil(numel(x)/20));   %Sample 1st 5% of points
        td = x(find(abs(sy - mean(y1)) > 3*std(y1),1));
            if isempty(td); td = 0; end
        sp = [c,a,b,td];     %Baseline, delay, set start
        sp = [sp(1), repmat(sp(2:end), 1, mn.n)];
        %Set coefficient bounds
        rx = range(x);	bl = [-ry, -ry, -ry, -rx];  
        bu = [ry,  ry,  ry,  rx];
        bl = [bl(1), repmat(bl(2:end), 1, mn.n)];
        bu = [bu(1), repmat(bu(2:end), 1, mn.n)];
        %Set both to identical for Two-term exponential ??
    case 'hill'     %Hill equation
        %   Coefs: c0, a1, ka1, h1
        c = sy(1);   a = sy(end) - c; h = 1;      %Baseline, amplitude, hill
        ind = find(sy > (sy(1) + sy(end))/2, 1);   %Index of half-max
        ka = x(ind);    sp = [c,a,ka,h];        %Set half-max and starts
        sp = [sp(1), repmat(sp(2:end), 1, mn.n)];
        %Set coefficient bounds
        rx = range(x);  bl = [min(y), -ry, -rx, 0  ]; 
        bu = [max(y),  ry,  rx, Inf]; 
        bl = [bl(1), repmat(bl(2:end), 1, mn.n)];
        bu = [bu(1), repmat(bu(2:end), 1, mn.n)];
        %Set both to identical for Two-term hill ??
    case 'gauss'
        %   Coefs: a1, b1, c1  Y = a1*exp(-((x-b1)/c1)^2)
        a1 = range(sy);  [mx,mi] = max(sy);   b1 = x(mi);   c1 = 1; %#ok<ASGLU>
        sp = [a1, b1, c1];  	sp = repmat(sp, 1, mn.n);
        %Set coefficient bounds
        bl = [-ry, min(x), 0  ];       bl = repmat(bl, 1, mn.n);
        bu = [ ry, max(x), Inf];       bu = repmat(bu, 1, mn.n);
    case 'fourier' 
        %   Coefs: a0, a1, b1, w  Y = a0+a1*cos(x*w)+b1*sin(x*w)
        a0 = mean(sy);  a1 = range(sy)/2; 	b1 = range(sy)/2;  
        w = 2*pi/range(x); 	sp = [a0, a1, b1, w]; 	
        sp = [sp(1), repmat(sp(2:end-1), 1, mn.n), sp(end)];
        %Set coefficient bounds
        bl = [min(y), -ry, -ry, 0];   
        bu = [max(y),  ry,  ry, 2*pi/min(abs(diff(x)))]; 
        bl = [bl(1), repmat(bl(2:end-1), 1, mn.n), bl(end)];        
        bu = [bu(1), repmat(bu(2:end-1), 1, mn.n), bu(end)];
end

%Concatenate outputs and return
mp = {bl, bu, sp};

end


