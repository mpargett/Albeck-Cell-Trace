%TF_PLOT
%   Plotting function for estimated transfer functions (provided complex as
%   a function of frequency)
%
%   [] = tf_plot(tf, fq, ch, x, y, t, opt)
%
%   Generates a plot for transfer functions, and input/output time series
%   data if provided.
%
%   


%Notes for new plotting function:
%   Basic plots: Gain, Phase, Coherence
%       plot w/ Confidence as on option
%   Call only Basics requested by param
%   Loop for overlays
%   Mean/Quantiles
%   

function out = ct_tfplot(atf, varargin)
%% Input parsing

%Default inputs (as a structure)
in = struct('fq', 0:0.5/(numel(atf)-1):0.5, ...
    'ch', [], 'x', [], 'y', [], 't', []);
opt = struct('xscale', 'log', 'plotgain', true, 'plotphase', true, ...
    'type', 'meanover', 'tu', 'au', 'fqu', 'Hz', 'lw', 1, 'cmap', 'o',...
    'cbs', []);

%Parse inputs to the input structure
inord = fieldnames(in);  nno = numel(inord);  nin = nargin-1;
for s = 1:min(nin,nno);   in.(inord{s}) = varargin{s};   end
%   Propagate options as needed
if nin > nno; ofn = fieldnames(varargin{nno+1});
    for s = 1:numel(ofn);  opt.(ofn{s}) = varargin{nno+1}.(ofn{s}); end
end

%Input matching for TF data
if isnumeric(atf);      atf = {atf};     end;  ntf = numel(atf);
nc = structfun(@isnumeric, in);  fnc = fieldnames(in);
for s = find(nc)';   in.(fnc{s}) = {in.(fnc{s})};    end
%   Ensure Row vector TFs etc.
atf = cellfun(@(x)reshape(x,1,numel(x)), atf, 'Un', 0);
in = structfun(@(z)cellfun(@(x)reshape(x,1,numel(x)),z,'Un',0),in,'Un',0);

%Input matching for frequency
if ntf > numel(in.fq);  in.fq(end+1:ntf) = in.fq(end);  end

%Check and set colormap as needed
if ischar(opt.cmap); opt.cmap = ct_colormaps(opt.cmap); end  %IF named color, get map
if isempty(opt.cmap); opt.cmap = colormap; close; end 
opt.cmap_x = linspace(0, 1, size(opt.cmap,1));	%Define colormap x-axis

%Generate plot flags from inputs
flg = {};   pio = false;
%   Regarding input/output
b = 0;
if ~all(cellfun(@isempty,in.x))
    b = b + 1; flg{b} = 'input';  pio = true; end
if ~all(cellfun(@isempty,in.y))
    b = b + 1; flg{b} = 'output'; pio = true; end
nio = b;    %Store number of IO panels
if isempty(in.t);  opt.tu = 'au';   end
%   Regarding TF
b = 0;
if opt.plotgain;     b = b + 1;  flg{nio+b} = 'gain';    end
if opt.plotphase;    b = b + 1;  flg{nio+b} = 'phase';   end
if ~all(cellfun(@isempty,in.ch));  b = b + 1;  flg{nio+b} = 'cohere';  end
npp = b;    %Store number of TF panels




%% Plot generation
fa = [];
if pio  %IF plotting inputs/outputs as well
    iof = figure(); set(iof, 'Position', [100 300 800 600]); fa = [fa,iof];
    for s = 1:nio
        switch flg{s}
            case 'input';    xx = in.t;     yy = in.x;
            case 'output';   xx = in.t;     yy = in.y;
        end
        ax = subplot(nio, 1, s);  opt.top = s == 1;  opt.bot = s == nio;
        sub_plotter(ax, flg{s}, xx, yy, opt);
    end
end

%Initialize figure
ff = figure();  set(ff, 'Position', [100 300 800 600]);
%Plots
for s = nio+(1:npp)  %For each plot to generate
    switch flg{s}
        case {'gain', 'phase'};  xx = in.fq;    yy = atf;
        case 'cohere';           figure(ff);   xx = in.fq;    yy = in.ch;
    end
    ax = subplot(npp, 1, s-nio); 
    opt.top = s-nio == 1;  opt.bot = s-nio == npp;
    yy = sub_plotter(ax, flg{s}, xx, yy, opt);
    switch flg{s}
        case 'gain';    out.gain = yy;  out.freq = xx;
        case 'phase';   out.phase = yy;
        case 'cohere';  out.coherence = yy;
    end
    %FIXME - show confidence bands (frequency) based on sample rate and
    %   time series length (overlay on gain and phase). What about
    %   coherence?  Show both bounds on coherence, or even on all?
    if ~isempty(opt.cbs);  sub_plotbnds(ax, flg{s}, opt.cbs);  end
end

out.fh = [fa, ff];

end


%% Subfunction: Plotter (flag for different plots)
function y = sub_plotter(ax, flg, x, y, p)

%Allow naive x-axis (initialize)
if isempty(x) || isempty(x{1})
    x = cellfun(@(yy)1:numel(yy),y,'UniformOutput', false); 
end

xscl = 'linear';  isfq = false;
ttl = 'Transfer function gain and phase';
switch lower(flg)
    case 'gain'
        y = cellfun(@(y)20.*log10(abs(y)), y, 'UniformOutput', false);
        ylbl = 'Gain (dB)';  xlbl = ['Freq. (',p.fqu,')'];
        xscl = p.xscale;   isfq = true;
    case 'phase'
        y = ct_getphase(y);     ylbl = 'Phase (deg.)';
        xscl = p.xscale;  xlbl = ['Freq. (',p.fqu,')'];  isfq = true;
    case 'cohere'; 	ylbl = 'Coherence';   xlbl = ['Freq. (',p.fqu,')'];
        xscl = p.xscale;    isfq = true;
    case 'input';  	ylbl = 'Input';     xlbl = ['Time (',p.tu,')'];
        ttl = 'Time Series Signals';
    case 'output';  ylbl = 'Output';    xlbl = ['Time (',p.tu,')'];
end

%Draw plots (overlays)
%   Normalize color intensities to use for ensemble plot
    ny = numel(y);  clr = linspace(0.5,0.8,ny)';  ndark = min(ny, 50);
    clr(1:ndark) = clr(1:ndark) + linspace(-0.1,0,ndark)';
%   Map color intensities to current colormap
    cvals = interp1(p.cmap_x, p.cmap, clr);
%   First, set up the axes and color order
set(ax, 'ColorOrder', cvals(end:-1:1,:));  hold(ax, 'on');

switch p.type
    case 'meanover'
        sub_plotens(ax, x, y, ny, p);  sub_plotmean(ax, x, y, p);
    case 'mean';        sub_plotmean(ax, x, y, p);
    case 'ensemble';    sub_plotens(ax, x, y, ny, p);
    case 'hist';        sub_plothist(ax, x, y, ny, p, isfq);
end

%Adjust axis and labeling
set(ax, 'xscale', xscl);       %Set x-axis scale ('log' or 'linear')
ylabel(ax, ylbl);                  %Set y-axis labeling
mnx = min( cellfun(@(z)min(z(z~=0)), x) );  mxx = max( cellfun(@max, x) );
set(ax,'XLim',[mnx, mxx]);

%If x-axis is frequency, add second axis for time per cycle
if isfq; box(ax, 'on');
    axy(1) = ax;  axy(2) = axes('Position', get(axy(1),'Position'),...
        'XAxisLocation', 'top', 'Color', 'none', 'YTick', []);
    xticks = 1./get(axy(1),'XTick');
%     xticks(1:end-1) = ceil(xticks(1:end-1));  
%     xticks(end) = floor(xticks(end));
    xtl = num2str(xticks(:),'%g');  
    set(axy(2), 'XLim', get(axy(1),'XLim'), 'XTick', 1./xticks, ...
        'XTickLabel', xtl, 'XScale', xscl);
    linkprop(axy, {'Position', 'XLim', 'XLimMode', 'YLim', 'YLimMode',...
        'XTick', 'XTickMode'});
    uistack(axy(1),'up');  %Move original axis to the top
else
end

%Apply Title and x-axis labeling
if p.top
    %Print title
    h_ttl = title(ttl);
elseif p.bot
    xlabel(ax, xlbl); %Apply x Label
end

end



%% Subfunction: Calculate mean and quantiles
function [mx, my] = sub_getmean(x, y)
%Interpolate for x and y, as needed
if iscell(x) && numel(x) > 1; 
    mnx = min(cellfun(@min, x));    mxx = max(cellfun(@max, x));
    lng = max(cellfun(@numel, x));
    mx = linspace(mnx, mxx, lng);
    y = cellfun(@(xi,yi)interp1(xi,yi,mx, 'linear', NaN), x, y, ...
        'UniformOutput', false);    y = cat(1, y{:});
else
    if iscell(y); y = y{1}; end
    if isempty(x); mx = 1:numel(y); elseif iscell(x); mx = x{1}; end
end

%Take mean
my = nanmean(y, 1);
%Get quantiles (0.25, 0.75)
my(2:3,:) = prctile(y, [25,75], 1);

end


%% Subfunctions for plotting
function [] = sub_plotens(ax, x, y, ny, p)
for s = 1:ny
    if isempty(x{s});   x{s} = 1:numel(y{s}); 	end
    plot(ax, x{s}, y{s}, 'LineWidth', p.lw);
end
end

function [] = sub_plotmean(ax, x, y, p)
    [mx, my] = sub_getmean(x, y);
    %Map dark values from colormap
    mqclr = interp1(p.cmap_x, p.cmap, 0.3);  mlow = mqclr < max(mqclr);
    mqclr = mqclr.*(mlow*0.5 + ~mlow*0.75);  %Darken and purify color
    %Plot
    plot(ax, mx, my(2:3,:), 'Color', mqclr, 'LineWidth', p.lw+1); %Quantiles
    mqclr = mqclr.*0.75; %Darken futher for the Mean plot
    plot(ax, mx, my(1,:),   'Color', mqclr, 'LineWidth', p.lw+2); %Mean
end

function [] = sub_plothist(ax, x, y, ny, p, isfq)
if isfq
    %Log scale and interpolate to same fqs
else
    %Assemble NaN padded array
    dx = min(x{1}(2:end) - x{1}(1:end-1));
    mxt = max(cellfun(@(xx)xx(end), x))/dx;
    yt = nan(ny, mxt);
    for s = 1:ny;  yt(s, x{s}/dx) = y{s}; end;    y = yt;
end
%Plot via CT_TRACKVIS
ct_trackvis(ax, 'histogram', y, 'cmap', p.cmap); 
end

function [] = sub_plotbnds(ax, flg, cbs)
hold(ax, 'on');
switch flg
    case 'gain';        pb = cbs(1,:);
    case 'phase';       pb = cbs(2,:);
    case 'cohere';      pb = [cbs(1,:),cbs(2,:)];
    otherwise; return
end
axy = axis(ax); yr = axy(3:4);
plot(ax, [pb;pb],  yr(:), 'r:');

end




