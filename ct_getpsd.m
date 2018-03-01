%CT_GETPSD
%   Calculate the Power Spectral Density of a signal
%
%   out = ct_getpsd(x, varargin)
%   
%   Returns a structure containing:
%       psd - Power Spectral Density estimates
%       f   - Frequency bins at which power is estimated
%       con - Confidence intervals (95%) for each frequency bin
%           Outputs match the input type (numeric for numeric,
%           multi-element structure for a cell array, structure with names
%           for a structure) 
%
%   Input options may be passed as field/value pairs
%       fq     - [1]    The sampling frequency of the signal
%       tunit  - ['s']  The time unit (matching sampling frequency)
%       iunit  - ['au'] The intensity unit (amplitude) of signal
%       meth   - ['adapt'] String indicating method to use:
%                'welch' for Welch's method, 
%                'adapt', 'eigen', 'unity' for Thomson multitaper method,
%                with adaptive, eigenvalue or unity weighting, respectively
%       center - [false]   Logical to center each signal (subtract mean)
%       mean   - [false]   Logical to produce a mean PSD
%       plot   - [false]   Logical to produce of plot of PSDs
%       plotlogp - [true]   Logical to plot Power in logspace
%       plotlogf - [true]   Logical to plot Frequency in logspace


function out = ct_getpsd(x, varargin)
%% INPUTS AND PROCESS PARAMETERS
%Default parameters
p.fq = 1;
p.tunit = 's';
p.iunit = 'au';
p.meth = 'welch';
p.center = false;
p.mean = false;
p.plot = false;
p.plotlogp = true;
p.plotlogf = true;


%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

%Check inputs and get number of cell trace samples
if isnumeric(x); nsig = 1; x = {x}; xc = 'numeric';
elseif iscell(x); nsig = numel(x); xc = 'cell';
elseif isstruct(x); xn = fieldnames(x); x = struct2cell(x); xc = 'struct';
end


%% PSD ESTIMATION
%Remove DC offset and 'center' each trace, if requested
if p.center; x = cellfun(@(xx)xx-mean(xx),x,'UniformOutput',false); end

psd = cell(nsig,1);  f = psd; con = psd; t = psd;
for s = 1:nsig
    nt = numel(x{s});              %Number of time points
    t{s} = (0:1:nt-1)./p.fq;  %Time domain indices
    if strcmpi(p.meth,'welch')  %IF Welch's method is preferred
        [psd{s}, f{s}, con{s}] = pwelch(x{s},[],[],[],p.fq);
    else
    %Get Transfer Function estimate (i.e. freq. response)
    [dpssE, dpssV] = dpss(nt, 4);
    [psd{s}, f{s}, con{s}] = pmtm(x{s}, dpssE, dpssV, ...
        [], p.fq, p.meth, 'ConfidenceLevel', 0.95); %#ok<PMTMCONF>
    end
end


%% PLOT
%   IF requested
if p.plot
    %Plot time domain signal(s)
    figure; subplot(2,1,1); hold on; title('Signal in Time Domain');
    for sp = 1:nsig; plot(t{sp},x{sp}, 'Color', [0,0.75,1]); end
    xlabel(['Time (',p.tunit,')']);
    ylabel(['Intensity (',p.iunit,')']);
    
    %Plot PSD(s)
    subplot(2,1,2); hold on; title('Power Spectral Density');  
    if p.plotlogp; p_psd = log10(abs(psd{sp})); p_con = log10(abs(con{sp}));
        yln = 'ln';
    else    p_psd = abs(psd{sp}); p_con = abs(con{sp}); yln = '';
    end
    for sp = 1:nsig
        plot(f{sp}, p_psd, 'Color', [0,0.75,1]);
        plot(f{sp}, p_psd, 'bo', 'MarkerSize', 2); 
        if nsig == 1; %IF only one, plot confidence intervals
        plot(f{sp}, p_con(:,1), 'Color', [0,0.75,0.75]);
        plot(f{sp}, p_con(:,2), 'Color', [0,0.75,0.75]);
        end
    end
    if p.plotlogf; set(gca,'xscale','log'); end
    axis tight;
    xlabel(['Frequency (1/',p.tunit,')']);
    ylabel(['Power (',yln,' Amp^2*',p.tunit,')']);
end


%% Match output style to input
switch xc
    case 'numeric'; out.psd = psd{1}; out.f = f{1}; out.con = con{1};
    case 'cell';    out.psd = psd;    out.f = f;    out.con = con;
    case 'struct'
        out.psd = cell2struct(psd, xn, 1);
        out.f   = cell2struct(f, xn, 1);
        out.con = cell2struct(con, xn, 1);
end

end