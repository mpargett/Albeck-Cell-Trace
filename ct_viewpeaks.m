%CT_VIEWPEAKS
%   Creates a plot to examine peaks as identified and described via
%   ct_pulseanalysis.
%
%   Usage: ct_viewpeaks(D, Z);
%          ct_viewpeaks(H, ...);
%               Generates plots on the figure or axis handle H
%           Additional inputs may be passed as parameter/value pairs.
%
%Inputs:
%   D - Data, as a numeric vector, matrix, or cell array
%   Z - The output structure from running ct_pulseanalysis on data D.
%   
%Parameters
%   nmax    - Maximum number of subplots, vertically
%   mmax    - Maximum number of subplots, horizontally
%   pos     - Position property of figure windows (e.g. to control size)

function [] = ct_viewpeaks(varargin)
%Default parameters
p.nmax = 5;     %Max N, number of axes vertically
p.mmax = 3;     %Max M, number of axes horizontally
p.pos = [100, 100, 1200, 800];
ah = [];    sp = true;   %Default empty axes handle and flag to subplot

%Manage passed handles
if ishandle(varargin{1}); ah = varargin{1};  varargin(1)  = []; end

%Collect inputs (data and pulse structure)
[d,z] = deal(varargin{1:2});  varargin(1:2) = [];
if ~iscell(d); d = num2cell(d, 2); end
nc = size(d,1); %Get size of input data

%Check validity if axes handle is provides
if ~isempty(ah); sp = false;
    assert(strcmpi(get(ah,'Type'),'axes'), 'Handle must be for axes.');
    assert(nc == 1, ['When providing an axes handle, only a single ',...
        'trace may be plotted.']);
end

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end


%% Plot each trace on a separate axis
nper = p.nmax.*p.mmax;  %Plots allowed per figure
for sc = 1:nc
    if sp
        %Get per figure axis index
        axi = mod(sc-1,nper) + 1;
        %Create new figure as needed
        if axi == 1; fh = figure(); set(fh, 'Position', p.pos); end
        %Create axes
        ah = subplot(p.nmax, p.mmax, axi);
    end
    
    %Plot data
    plot(ah, d{sc}, 'Color', [0, 0.8, 0.6]); hold on;
    
    %Indicate peaks
    ax = axis;  yr = diff(ax(3:4));
    plot(ah, z(sc).pkpos, d{sc}(z(sc).pkpos) + yr*0.025,'kv');
    
    %Draw peak features
    for s = 1:numel(z(sc).pkpos)
        %Show peak amplitudes (as a line from peak to base)
        plot(ah, z(sc).mpos([s,s]), z(sc).mean(s)-[0,z(sc).amp_peak(s)], 'b-');
        %Show peak durations (as horizontal line at half-height)
        plot(ah, z(sc).mpos(s)+z(sc).dur(s).*[-0.5,0.5], ...
            z(sc).mean([s,s]) - z(sc).amp_mean(s)./10, 'b-');
        %Show rise and fall regions
        rts = z(sc).mpos(s) - 0.5*z(sc).dur(s) +[-z(sc).rise(s), 0];  %Rise time points
        fts = z(sc).mpos(s) + 0.5*z(sc).dur(s) +[0, z(sc).fall(s)];   %Fall time points
        plot(ah, rts, d{sc}(ceil(rts)), 'r:');  
        plot(fts, d{sc}(floor(fts)), 'r:');
    end
end


end