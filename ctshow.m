%CTSHOW
%   Cell Trace SHOW for image data.  Trims data for expected usable
%   segments and plots heatmaps and overlayed traces for visualization.
%
%   dout = ctshow(in)
%   where in is a structure containing valcube as a field.
%   dout provides a structure containing the 'valid' data ranges, time
%   index and Ratio, RFP data used (included interpolation across small
%   gaps)
%
%       e.g. usage
%   in = load(path_to_processed_data_file);
%   ctshow(in);

function [dout] = ctshow(in,pin)

%Plot options
lw = 2;

%Default parameters
%   Get channel indices from provided order, if available
defc = {'EKAR','RFP'};  defm = {'cfp.?cyt[^/]*/yfp.?cyt', '^rfp.?nuc[^/]*'};
if isfield(in, 'vcorder')
    for s = 1:numel(defc)
        p.chan{s} = defc{s};
        p.ind{s} = find(~cellfun('isempty', regexpi(in.vcorder, defm{s})));
        if numel(p.ind{s}) ~= 1; p.ind(s) = {}; continue; end
    end
else    %Otherwise, use old defaults
    p.ind = [9,3];                              %Channel indices (valcube)
    p.chan = {'EKAR','RFP'};                    %Channel names
end
p.fun = {@(x)x, @(x)x};  %Input/Output functions (as is)
p.fretfactor = 1;
p.gapmax = 2; p.dlengthmin = 100;           %Trimming parameters
p.tsamp = 7.5;                              %Time step, minutes
p.startonly = true;

%Apply input parameters
if exist('pin','var')
pn = fieldnames(pin); for s = 1:numel(pn); p.(pn{s}) = pin.(pn{s}); end
end

%Check on dlengthmin (minimum value = 2)
if p.dlengthmin < 2
    p.dlengthmin = 2; 
    warning('Minimum value for dlengthmin is 2. Using dlengthmin = 2');
end

nchan = length(p.ind);

%Parse input function names, if provided
yt = cell(nchan,1);
for sn = 1:nchan;   yt{sn} = 'Intensity';  %#ok<ALIGN>
    if ischar(p.fun{sn});
        switch lower(p.fun{sn}) %#ok<ALIGN>
            case 'fret';    p.fun{sn} = @(x,y)1 - x./y;
            case 'local';   p.fun{sn} = @(x,y)x./y;
            case 'invert';  p.fun{sn} = @(x)bsxfun(@minus,2*mean(x,2),x);
            case 'fret_est'
                p.fun{sn} = @(x)1 - x./p.fretfactor;
                yt{sn} = 'Est. FRET Ratio';
end; end; end

%Collect and pretreat data (as needed)
if isfield(in,'data');
    d = in;
elseif isfield(in, 'valcube')
    for s = 1:nchan
        %Get desired data, checking for multiple inputs
        if iscell(p.ind); nins = numel(p.ind{s}); din = {};
        for ss = 1:nins
            din{ss} = in.valcube(:,:,p.ind{s}(ss)); %#ok<AGROW>
        end
        else din{1} = in.valcube(:,:,p.ind(s));
        end
        %Convert zeros to NaNs for robust recognition
        din = cellfun(@(x)nanify(x), din, 'UniformOutput', false);
        %Process with data functions
        d.(p.chan{s}) = p.fun{s}(din{:});
    end
    
    %Trim data
    trimp = struct('gapmax',p.gapmax,'dlengthmin',p.dlengthmin,...
        'startonly',p.startonly);
    d = ct_trimdata(d, trimp);
end
dout = d;

%% Plot
for s = 1:nchan
    %Assemble full maxtrix version of data
    ndat = length(d.data); 
    ntime = d.time(:,2) - d.time(:,1) + 1;
    ntmax = max(d.time(:,2));  dtemp = zeros(ndat,ntmax);
    for ss = 1:ndat
        dtemp(ss, d.time(ss,1):d.time(ss,2) ) = d.data(ss).(p.chan{s});
    end
    
    %Plot heatmap of all cell traces
    figure();  subplot(2,1,1);    
    imagesc(dtemp);  colormap(gca,hot); 
    %   Get Data bounds
    mVal = min([d.data.(p.chan{s})]);
    xVal = max([d.data.(p.chan{s})]);
    %   Set color axis for best visibility
    caxis([mVal, xVal]);
    ps = get(gca,'Position'); 
    set(gca,'Position', ps + [0,-0.05,0,0.05],...
        'XTickLabel', [], 'YTickLabel', [], 'XTick', [], 'YTick', []);
    ylabel('Single Cells');  title(p.chan{s}); box on;
    
    %Plot ensemble cell traces
    subplot(2,1,2); hold on;
%     ntdev = zeros(ndat,1);
%     for ss = 1:ndat
%     %   Determine coloration net deviation
%     ntdev(ss) = sqrt( mean((d.data(ss).(p.chan{s})(:,2:end) ...
%                         - d.data(ss).(p.chan{s})(:,1:end-1)).^2) )...
%                 ./(mean(d.data(ss).(p.chan{s}),2)+1);
%     end
    % Normalize colors
%     clr = ntdev - min(ntdev); clr = 0.9 - 0.6*clr./(max(clr) + eps);
%     [clrsort, plotorder] = sort(clr,1,'descend');
    %HERE, RANDOM COLORATION (grayscale)
    clr = rand(ndat,1);  clr = 0.9 - 0.6*clr./(max(clr) + eps);
    %   Now Plot
    for ss = 1:ndat
        xx = ( d.time(ss,1):d.time(ss,2) )*p.tsamp;
        plot(xx, d.data(ss).(p.chan{s}), 'Color', ...
            [clr(ss),clr(ss),clr(ss)], 'LineWidth',lw);
    end
    
    %Plot mean value at each time point
    %   Get mean value
    mv = sum(dtemp,1)./sum(dtemp~=0,1); 
    %   Get quantiles (25% and 75%)
    qt = quantile(nanify(dtemp), [0.25, 0.75], 1);
    %   Add to plot
    xx = (0:ntmax-1)*p.tsamp;
    plot( xx, qt, 'Color', [1,1,0], 'LineWidth', lw-1);
    plot( xx, mv, 'Color', [1,0.9,0], 'LineWidth', lw);
    
    %Set axis for best visibility
    mt = max(ntime*p.tsamp);
    axis([-0.002*mt, mt, mVal, xVal]);
    ps = get(gca,'Position'); set(gca,'Position', ps + [0,0,0,0.05]);
    xlabel('Time (min)'); ylabel(yt{s});
    set(gcf,'Position',[100, 350, 800, 600]);
    
    %Store mean values
    dout.mean.(p.chan{s}) = mv;
end

end


%% NaNify
function x = nanify(x)

if isnumeric(x)
    x(x == 0) = NaN;
elseif iscell(x)        %Recursive for cells
    x = cellfun(@(z)nanify(z), x, 'UniformOutput', false);
elseif isstruct(x)      %Recursive for structures
    x = structfun(@(z)nanify(z), x, 'UniformOutput', false);
end

end