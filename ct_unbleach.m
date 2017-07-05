%CT_UNBLEACH
%	Remove estimated photobleaching artifact from Cell Traces
%	Uses a single exponential fit to estimate bleaching rate

function [ct, pv] = ct_unbleach(ct, pin)

%Establish parameters
%   Default settings
p.ix = 13;  p.iy = 14;
p.ic.CFP = [1,4]; p.ic.YFP = [2,5];  p.ic.RFP = [3,6];
p.tsamp = 1;  p.tir = [1,size(ct,2)];
%   Apply input parameters
if exist('pin','var') && isstruct(pin) && ~isempty(pin)  
    pn = fieldnames(pin);   for s = pn(:)'; p.(s{1}) = pin.(s{1}); end
end

imn = [cell2mat( struct2cell(p.ic)' ), p.ix, p.iy];
cf = fieldnames(p.ic);    nc = numel(cf);

%% Get mean value over traces
%	Ensure only tracks (and times) with valid entries and coords are used
ct_valid = all(~isnan(ct(:,:,imn)) & ct(:,:,imn) ~= 0,3);
%   Retain only traces that exist at the start of the time range
ct_valid(~ct_valid(:,p.tir(1)),:,:) = false;  
%   Eliminate empty rows
userow = any(ct_valid, 2);  nr = nnz(userow);
%   Elimnate out of time range elements
ctt = ct(userow,p.tir(1):p.tir(2),:);  
ct_valid = ct_valid(userow,p.tir(1):p.tir(2));
%   Negate remaining invalid entries
ctt(repmat(~ct_valid, [1,1,size(ct,3)])) = 0;
%   Get time range to consider from remaining tracks
nt = find(any(ct_valid,1),1,'last');
tvec = p.tsamp.*( 0:(nt-1) );  %-1 to make first value 0

%Combine populations from the same channel
ctt = num2cell(ctt, [1,2]);
ctt = structfun(@(x)cat(1, ctt{x}), p.ic, 'UniformOutput', false);
%Normalize to first image intensity
ctt = structfun( @(x)bsxfun( @rdivide, x, x(:,1) ), ctt, ...
    'UniformOutput', false);
%Evaluate means
nctv = sum(ct_valid,1)./nr;
mv = structfun( @(x) bsxfun( @rdivide, sum(x,1), nctv*size(x,1) ), ctt, ...
    'UniformOutput', false); 

%Eliminate gaps in the data by interpolation, in any are present
%   Interpolate any gaps in data
gaps = structfun( @isnan, mv, 'UniformOutput', false);
gaps = any(cell2mat(struct2cell(gaps)),1);
if any(gaps(:))
    iv = (1:nt)';  %Index vector over time
    mvt = cell2mat(struct2cell(mv));
    mvt(:,gaps) = interp1(iv(~gaps), shiftdim(mvt(:,~gaps),1), iv(gaps), ...
        'linear', 'extrap')';
    mv = cell2struct( num2cell(mvt, 2), cf, 1);  clear mvt;
end

    
%% Fit mean values to single exponential decay models
%   Define model
ff = @(a,k,t)a + (1-a).*exp(-exp(k).*t);
%   Define parameter bounds
LBs = [0, -20];       UBs = [1, 1];
%Fit by minimized sqared Euclidean error
y = nan(nc,2);     rsq = zeros(nc,1);
fmcopt = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
    'Display', 'off');
for s = 1:nc
    [y(s,:)] = fmincon(...
        @(x)sum((ff(x(1),x(2),tvec) - mv.(cf{s})).^2) + nt*8e-4*(x(1)).^2, ...
        [0,-5], [],[], [],[], LBs, UBs, [], fmcopt);
    %   For Refernce, the L2 regularization using a factor of 8e-4
    %   	estimates 1/100 the cost of a poor fit (i.e. regularization
    %   	would compete with a good fit, not a poor one.
    %Get validty metrics (correlation coefficient here)
    rsq(s) = corr(shiftdim(mv.(cf{s}),1), ...
                  shiftdim( ff(y(s,1),y(s,2),tvec), 1)).^2;
end

usemodel = rsq > 0.5 & y(:,1) < 0.99 & exp(y(:,2))*tvec(end) > 0.01;
%Check for validity of model


%% Perform correction on data
%   Reset time vector to encompass full data range
tvec = p.tsamp.*(0:(size(ct,2)-1));  %-1 to get time value from index
pv.Mod = ff;     pv.tsamp = p.tsamp;
for s = 1:nc
    if usemodel(s)
    for ss = 1:length(p.ic.(cf{s}))
        ct(:,:,p.ic.(cf{s})(ss)) = bsxfun( @rdivide, ...
            ct(:,:,p.ic.(cf{s})(ss)), ff(y(s,1), y(s,2), tvec) );
    end
    end
    %Provide exponential models and fitness, as desired
    pv.(cf{s}) = struct('a',y(s,1),'k',y(s,2), ...
                        'rsq',rsq(s), 'used',usemodel(s));
end



end