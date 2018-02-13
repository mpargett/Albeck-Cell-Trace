%CT_INFERLTI
%   Perform inference of LTI systems for Celltracer time series data.
%   Constructs Transfer Function (TF) representations in frequency space,
%   from pairs of time series data, or infers TF parameters (zeros,poles)
%   from TFs, comparing with time series data to calculate error.
%
%   [out] = ct_inferlti('est', x, y, 'Param', Value, ...)
%       estimates TFs from time series data (x and y).  x and y must be
%       paired cell arrays of data vectors.  Additional parameters are:
%       Center - [log] TRUE to center all data (subtract mean)
%       Meth - [char] Method of estimation to use ('adapt',
%           'unity', 'eigen' or 'welch' (default)).  See pmtm.
%       fqs - [num] Sampling frequency (i.e. 1/t_samp)
%       MeanTF - [log] TRUE to calculate a mean TF from CPSDs
%
%   [out] = ct_inferlti('fit', sys, x, y, 'Param', Value, ...)
%       infers parameters for the TFs supplied in sys, the output structure
%       from ct_inferlti('est',...), considering various numbers of
%       parameters and choosing the best parsimonius fit by AIC.
%       Additional parameters are:
%       nbmax - [num] Max number of Zeros to consider
%       namax - [num] Max number of Poles to consider
%

%FIXME: Calculate time domain residuals - toward unknown inputs etc.
%FIXME: Calculate reasonable confidence bounds on the frequency range in
%   which estimates may be unbiased by windowing/tapers etc.  Resort to
%   empirical derivation from simulated data if needed.

function out = ct_inferlti(flag, varargin)
nin = length(varargin);   	%Get number of inputs

%% Perform major LTI inference operations
switch lower(flag)
    case 'est'  %Estimate LTI Transfer Function for each Time Series (x,y)
        %Parse inputs
        [x, y] = deal(varargin{1:2});
        %   Decalre default parameters
        p.center = false;   p.meth = 'welch';   
        p.meantf = false;   p.fqs = 1;
        %   Split remaining pairs to a structure
        for s = 3:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end
        %Call LTI estimation routine
        out = sub_inferLTI(x, y, p);
    case 'fit'  %Fit LTIs in a constrained parameter space
        %Parse inputs
        [sys, x, y] = deal(varargin{1:3});
        %   Decalre default parameters
        p.nbmax     = 0;   	%Maximum TF numerator order   (~# of zeros)
        p.namax     = 10;	%Maximum TF denominator order (~# of poles)
        p.fqs = 1;
        %   Split remaining pairs to a structure
        for s = 4:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end
        %Call LTI parameter inference routine
        out = sub_fitLTIpar(sys, x, y, p);
end

end


%% Unconstrained LTI inference
function q = sub_inferLTI(x, y, p)
%Check input format and get number of cell trace samples
if isnumeric(x);    nsig = 1;  x = {x};  y = {y}; 
else                nsig = numel(x); 
    if isstruct(x);  x = struct2cell(x);  y = struct2cell(y);   end
end

%Remove DC offset and 'center' each trace, if requested
if p.center; x = cellfun(@(xx)xx-mean(xx),x,'UniformOutput',false); end
 
% Get TF Estimates for each trace
%   Use Multi-Taper averaging method, expecting non-periodic signals and
%   relatively short duration (multiple windows to be different signals,
%   not just different noise) 
[fq, tf_e, cn, ch] = deal(cell(nsig,1));  %Initialize
for s = 1:nsig  %FOR each signal
    nt = numel(x{s});              %Number of time points
    %   Check for size consistency (skip if bad)
    if nt ~= numel(y{s}); warning(['Input vectors (x, y) ',...
            'must be equal length.']);  continue;  end
    %Perform TF Estimation, per method specified
    if strcmpi(p.meth, 'welch')
        %   Welch's method based versions
        [tf_e{s}, fq{s}] = tfestimate(x{s}, y{s}, [],[],[], p.fqs);
        [ch{s}] = mscohere(x{s}, y{s}, [],[],[], p.fqs);
%         %Backup procedures, manually performing Welch
%         [psd_mt(s).cross, f_mt{s}] = cpsd(x{s},y{s},[],[],[],fq_samp);
%         [psd_mt(s).auto1] = pwelch(x{s},[],[],[],fq_samp);
%         [psd_mt(s).auto2] = pwelch(y{s},[],[],[],fq_samp);
    elseif strcmpi(p.meth, 'exact')
        %Get exact (overfit) TF
        nfft = numel(x{s});  nfft = nfft + rem(nfft,2);
%         nfft = 2^nextpow2(numel(x{s}));
        tf_e{s} = fft(y{s},nfft)./fft(x{s},nfft);
        tf_e{s} = tf_e{s}(1:(nfft/2 + 1))*2;
        fq{s} = p.fqs*(0:(nfft/2))/nfft;   ch{s} = [];   cn{s} = [];
    else
        %Get Transfer Function estimate (i.e. freq. response)
        [dpssE, dpssV] = dpss(nt, 4);  %Get orthogonal bases (Slepian)
        %   Calculate PSDs from MultiTaper base
        [psd_e(s), con_e(s), fq{s}] = pmtm_mod([x(s),y(s)], dpssE, ...
            dpssV, [], p.fqs, p.meth); %#ok<AGROW>
        %   Calculate estimated TF
        tf_e{s} = conj(psd_e(s).cross)./psd_e(s).auto1;
        %   Confidence estimates
        cn{s} = conj(con_e(s).cross)./con_e(s).auto1(:,end:-1:1);
        %   Get Coherence estimate (~confidence)
        ch{s} = (abs(psd_e(s).cross).^2)./(psd_e(s).auto1.*psd_e(s).auto2);
    end
end

%Pack TF values into output structure (ordered as were input signals)
q = cell2struct([fq, tf_e, cn, ch], {'fq', 'tf', 'cn', 'ch'}, 2);

%Mean Transfer Function
if p.meantf && ~strcmpi(p.meth, 'exact')
    %Define uniform frequency coordinates
    nf = min(cellfun('length', fq));  fmin = max(cellfun(@(x)x(1),fq));
    fmax = min((cellfun(@(x)x(end),fq))); f_mn = linspace(fmin,fmax,nf);
    %Interpolate PSDs to same frequencies (expect little difference)
    mn_cross = cellfun(@(x,y)interp1(x,y,f_mn), fq, {psd_e.cross}', ...
        'UniformOutput', false);  mn_cross = mean(cat(3, mn_cross{:}),3);
    mn_auto1 = cellfun(@(x,y)interp1(x,y,f_mn), fq, {psd_e.auto1}', ...
        'UniformOutput', false);  mn_auto1 = mean(cat(3, mn_auto1{:}),3);
    mn_auto2 = cellfun(@(x,y)interp1(x,y,f_mn), fq, {psd_e.auto2}', ...
        'UniformOutput', false);  mn_auto2 = mean(cat(3, mn_auto2{:}),3);
    %Calculate mean Transfer Function and Coherence
    tf_mean = conj(mn_cross)./mn_auto1;
    ch_mean = (abs(mn_cross).^2)./(mn_auto1.*mn_auto2);
    %   Pack into output structure
    q(nsig+1).fq  = f_mn;        q(nsig+1).tf = tf_mean;    
    q(nsig+1).ch = ch_mean;
end

end


%% Parameterized LTI fitting
function lfit = sub_fitLTIpar(sys, x, y, p)
%LTI fitting
%Fit parsimonious form for each TF (parsimonious over all of them?)
%   (Perform a search for coefficients, or zeros/poles, that minimize the
%       prediction error.  This is in accordance with methods used in the
%       System Identification Toolbox.)
%   Use BODE(TF(parameters),Frequncies), or FREQRESP(TF,Freq) to deliver
%   outputs and get error between predicted and measured response.
%DO: Loop TF ID with increasing order, until deltaAIC values based on the
%prediction errors fail to meet the threshold.

if isnumeric(x); x = {x}; end;   if isnumeric(y); y = {y}; end;
%Get number of signals
nsig = numel(sys);
%Get number of models to try per signal
nm = sum(1 + p.namax - (0:p.nbmax)) - 1;    

warning('Off', 'Control:analysis:LsimStartTime');
lfit = struct('tf',cell(nsig,1),'akaikewt',cell(nsig,1)); tf_p = cell(1,2);
for s = 1:nsig  %FOR each signal
    %Initialize per-signal variables
    SSE = zeros(nm, 1);  TFtemp = cell(nm,2); np = zeros(nm,1); ind = 0;
    for sb = 0:p.nbmax       %FOR different numbers of zeros allowed
        for sa = sb:p.namax  %FOR different numbers of poles allowed
            ind = ind + 1;  %Increment storage index
            %Fit TF equation to TF estimate
            
            fq_wts = [];
            %FIXME. IF desire frequency weights, apply here.
            %Define fitting weights
            %         fq_wts = (1:-0.9/(length(sys(s).f)-1):0.1)';
            %         fq_wts = exp(5*fq_wts)/exp(5); %Downweight high freqs
            
            [TFtemp{ind,1:2}] = invfreqz(sys(s).tf, ...
                sys(s).fq*pi/(p.fqs/2), sb, sa, fq_wts, 1e3);
            %       Uses normalized frequencies [0, pi]
            %Simulate time response of inferred system
            nt = numel(x{s});              %Number of time points
            t_samp = (0:1:nt-1)./p.fqs;  %Time domain indices
            tresp = lsim(  tf(TFtemp{ind,1:2}, 1/p.fqs),...
                [ones(1,nt-1)*mean(x{s}(1:10)), x{s}]', ...
                [-t_samp(end:-1:2),t_samp]  );
            %       Includes 'burn-in' period to rise to initial condition
            
            %Calculate AIC and Weights
            SSE(ind) = sum( ( tresp(nt:end) - y{s}' ).^2 );  %Sum(err^2)
            np(ind) = sb+sa;    %Number of parameters for this iteration
            %FIXME: Move AIC out of loop (or have termination criteria)
            %   Calculate new AICc values (all trials thus far)
            AICc = nt.*log(SSE(1:(ind))./nt) ...
                + 2*np(1:ind).*(np(1:ind) + 1)./(nt - np(1:ind) - 1);
            dAIC = AICc - min(AICc);                %delta AICc
            Awt = exp(-dAIC/2)./sum(exp(-dAIC/2));  %Akaike weights
        end
    end
    [mx, mi] = max(Awt);                %Find Max Akaike Wt
    [tf_p{1,1:2}] = TFtemp{mi,1:2};   %Keep TF from max Wt
    
    %Fill output structure
    lfit(s).tf = cell2struct(tf_p, {'b','a'}, 2);
    lfit(s).akaikewt = Awt;
    lfit(s).na = p.nbmax:sa;    lfit(s).nb = p.nbmax;  %???
    lfit(s).sse = SSE(mi);
end

%Akaike Information Criterion definitions:
%   AIC = ns.*log(SSE/ns) + 2*np;  %ns: # of samples; np: # of parameters
%   AICc = AIC + 2*np*(np + 1)/(ns - np - 1);  %Corrects for sample size
%   dAIC = AIC - minAIC;
%   AICwt = exp(-dAIC/2)/sum(exp(-dAIC/2));
%   Evidence ratio = AICwt1/AICwt2;  %1 is ... more likely than 2
%       Per Burnham KP, Anderson DR. Model selection and multimodel
%       inference: a practical information-theoretic approach. Springer,
%       2002.

%FIXME.  How to deal with different order fits for different cells???

%Get mean TF and Distribution Statistics
%   i.e. fit distributions of TF coefficients (Gaussian, perhaps) and
%   provide mean and variance values
end



