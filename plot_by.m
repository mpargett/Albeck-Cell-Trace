% PLOTBY
% This code plots imaging experiment data by celltype or by treatment type.
%
% INPUTS
%   plotby - either 'celltype' or 'treatment'. If set to 'celltype' function
%            will create a subplot for each cell type entered in data sheet
%            and plot all treatments on each subplot. 'treatment' creates a
%            subplot for each treatment.
%
%   data - should be a cell array of XYs containg .data.channel. This is
%          the output of the dataproc function.
%
%   channel - name of the channel you wish to plot as a string. Ex: 'ekar'
%
%   pmd - this is the first output of iman_readdatasheet. It is a structure
%         containging plate map data from your excel sheet.
% OPTIONAL INPUTS
%   lines - defualt is off. Set to true to plot a vertica line at treatment
%   times. Tx times must be input into data sheet as @##tp 
%
%   subset - only plot a subset of treatments or celltypes. For example, if
%   you want to plot only wells that were treated with EGF, subset = 'EGF'.
%
%
% EXAMPLE
% plot_by('celltype',d,'ekar',pmd,'lines',1)
%
% -------------------------------------------------------------------------
% VERSION 2.0
function [cellfn, txs] = plot_by(plotby,data,channel,pmd,varargin)
%% Check input
plotby = lower(plotby);

if any(strcmp(plotby,{'treatment','treatments','tx','txs'}))
    plotby = 'treatment';
elseif any(strcmp(plotby,{'cell','celltype','celltypes'}))
    plotby = 'celltype';
else error('Accepted plotby types: treatment or celltype. plotby must be a string')
end
%% Input parsing
p.lines = []; p.subset = []; p.tsamp = 1; p.hours = 0; p.dash = 1;
%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

if ~isempty(p.subset) && ~iscell(p.subset)
    p.subset = {p.subset};
end
%% get cell type names
c1 = pmd.Cell(:,:,1);
c2 = c1(cellfun(@isstr,c1)); celltypes = unique(c2); % find unique names of celltypes
celltypes = celltypes(~cellfun(@isempty,celltypes)); % discard empty name fields
% get rid of @ density if present for fieldname
cellfn = cellfun(@(x)x{1},regexp(celltypes,'@','split'),'un',0);
% remove any spaces in name
cellfn = arrayfun(@(x)regexprep(cellfn{x},'[^a-zA-Z0-9]',''),1:size(cellfn),'un',0);
% make list of xys for each celltype
for s = 1:numel(cellfn)
    %   Assemble list of xys with the current cname
    idx.(cellfn{s}) = [pmd.xy{strcmp(celltypes{s},c1)}];
end

%% Find good xys that actually have data
goodxy = false(1,max([pmd.xy{:}])); 
goodxy(~cellfun(@isempty,data)) = deal([true]);
%% get treatments and find unique combos
catTxnames = {'pTx','Tx1','Tx2','Tx3'}; % possible treatment fields
[~,I]= intersect(catTxnames, fieldnames(pmd)); % tx fields in experiment
I = sort(I); catTxnames = catTxnames(I);

Txcat = cell(size(pmd.Cell,1),size(pmd.Cell,2)); % initialize cell array
linetp = cell(size(pmd.Cell,1),size(pmd.Cell,2)); % initialize cell array

for s = 1:numel(catTxnames) % cat each treatment in a loop
    Txcat = cat(3,Txcat, cellfun(@(x)cat(2,x{:}),...
        cat(3,num2cell(pmd.(catTxnames{s})(:,:,1:3),3)),'Un',0));
    linetp = cat(3,linetp, cellfun(@(x)cat(2,x{:}),...
        cat(3,num2cell(pmd.(catTxnames{s})(:,:,4),3)),'Un',0));
end
linetp = linetp(:,:,2:end); % remove empty layer from initialization
linetp = num2cell(linetp,3);
% color sequence for treatment lines
linecolor = parula(numel(catTxnames)+1);

for s=1:numel(Txcat)
   if strcmp(Txcat{s},'')
       Txcat{s} = NaN;
   end
end

stridx = cellfun(@isstr,Txcat);
spacer =  repmat({'+'},size(Txcat)); % add + between txs for readability
spacer(~stridx) = {[]};

% add + between treatment names
Txcat2 = Txcat(:,:,2);
for s = 2:size(Txcat,3)-1
    Txcat2 = cat(3,Txcat2,spacer(:,:,s+1),Txcat(:,:,s+1));
end
% cat Tx name and space together into one long string
Txcat2 = cellfun(@(x)cat(2,x{:}),cat(3,num2cell(Txcat2,3)),'Un',0);
stridx = cellfun(@isstr,Txcat2);
treatments = unique(Txcat2(cellfun(@isstr,Txcat2)));

% replace spaces and / so txs can be used as fieldnames
for s = 1:numel(Txcat2)
    if isstr(Txcat2{s})
    Txcat2{s} = regexprep(Txcat2{s},char(0),'');
    Txcat2{s} = regexprep(Txcat2{s},'(^|\s|+)','');
    Txcat2{s} = regexprep(Txcat2{s},'[^a-zA-Z0-9]','_');
    end
end
% find the unique treatment combos
txs = unique(Txcat2(cellfun(@isstr,Txcat2)));
% restrict to subset
if ~isempty(p.subset)
    txI = ~cellfun(@isempty,(regexp(txs,strjoin(p.subset,'|'))));
    ctI = ~cellfun(@isempty,(regexp(cellfn,strjoin(p.subset,'|'))));
   if any(txI);
       txs = txs(txI);
       treatments = treatments(txI);
   end
   if any(ctI)
       cellfn = cellfn(ctI);
   end
    
end

% beautify treatment names for labels
% remove extra spaces in name caused by NaNs to make pretty names
treatments = arrayfun(@(x)regexprep(treatments{x},char(0),''),1:size(treatments),'un',0)';
treatments = arrayfun(@(x)regexprep(treatments{x},'(^|\s)',''),1:size(treatments),'un',0)';
% remove double +s
treatments = arrayfun(@(x)regexprep(treatments{x},'++','+'),1:size(treatments),'un',0)';
% replace _ with ' ' 
treatments = arrayfun(@(x)regexprep(treatments{x},'_',' '),1:size(treatments),'un',0)';


% vector of xys for packing txs into idx structure.

xymat = pmd.xy(stridx);
% now pack concatinated txs into the idx format 
for s = 1:numel(txs)
idx.(txs{s}) = cat(2,xymat{strcmp(txs{s},Txcat2(stridx))}); % find xy
% idx.(txs{s}) = 
end

linetp = linetp(stridx);


%% Plotting
switch plotby
    % Plot by Cell type
    case 'celltype'
        titlevec = cellfn;
        cmap = rainbow(numel(txs));
        legname = treatments;
        [TXX,nrow,ncol] = main_plotting_func(data,channel,cellfn,txs,idx,...
            xymat,cmap,linetp,titlevec,legname,goodxy,p);
    case 'treatment'
        titlevec = treatments; % use beautified tx names
        cmap = rainbow(numel(cellfn)+1);
        legname = cellfn;
        [TXX,nrow,ncol] = main_plotting_func(data,channel,txs,cellfn,idx,...
            xymat,cmap,linetp,titlevec,legname,goodxy,p);
        
end

% if requested, plot treatment lines
if p.lines
    for s = 1:numel(TXX)
        subplot(nrow,ncol,s),
        lnplot(TXX{s})
    end
end
end
%% Main plotting function 
function [TXX,nrow,ncol] = main_plotting_func(data,channel,segv1,segv2,...
    idx,xymat,cmap,linetp,titlevec,legname,goodxy,p)

  [nrow,ncol,pos] = subplot_sizer(numel(segv1)); % optimal subplot dimensions
  
        figure,
        for s = 1:numel(segv1)
            subplot(nrow,ncol,s), hold on
            sc = 0; legvec = [];
            % make dashed lines if there are a lot of lines
            if numel(segv2) > 5 && p.dash
                linestyle(1:2:numel(segv2)) = deal({'-.'});
                linestyle(2:2:numel(segv2)) = deal({'-'});
            else
                linestyle(1:numel(segv2)) = deal({'-'});
            end
            
            for sa = 1:numel(segv2) % plot line for each unique treatment/cell
                xy = intersect(idx.(segv1{s}),idx.(segv2{sa})); % find xy for treatment/celltype
                % loop through xy
                for sb = 1:length(xy) % plot line for each xy
                    if goodxy(xy(sb))% make sure there is data to plot.
                    % store tx times for each xy for plotting lines later
                 
                    sc = sc+1; % index for saving individual line objects for input to legend.
                    if sb == 1
                        legvec = [legvec,1]; % identify unique legend entries.
                        txx{sa} = str2double(linetp{...
                            cell2mat(cellfun(@(x)any(x == xy(sb)),xymat,'un',0))});
                    else
                        legvec = [legvec,0];
                    end
                    h{sc} = plot(nanmean(data{xy(sb)}.data.(channel)),...
                        'Color',cmap(sa,:),'LineStyle',linestyle{sa},'LineWidth',1.1 );
                    else
                        warning(['XY ',num2str(xy(sb)),' contains no data.'])
                    end
                end
            end
            
            % store tx times for plotting lines later
            TXX{s} = unique(cell2mat(txx));
            TXX{s} = TXX{s}(~isnan(TXX{s})); %#ok<*AGROW>
            clear txx
            
            % axis and title
            axis tight
            title(titlevec{s})
            
        end
        
        legvec = logical(legvec);
        subplot_standardizer(gcf) % make axis limits the same for all plots
        set(gcf,'Position',pos) % figure bigger
        extraleg(legname,[h{legvec}],nrow,ncol,s) % make legend in extra subplot
end
        
%% function to select figure dimensions
function [nrow,ncol,pos] = subplot_sizer(numplots)
% determines optimal subplot sizing
numplots = numplots+1; % extra plot for giant legend
if numplots > 6
    ncol = ceil(sqrt(numplots));
    nrow = ceil(numplots/ncol);
else
    ncol = numplots; nrow = 1;
end
hp = ncol*450; if hp > 1920, hp = 1920; end
vp = nrow * 250; if vp > 964, vp = 964; end
pos = [0,0,hp,vp];
end
%% function to set same x and y axis limits for all plots
function [] = subplot_standardizer(fig)

sp = get(fig, 'Children');

for i = 1:numel(sp)
    xlim_min(i) = sp(i).XLim(1);
    xlim_max(i) = sp(i).XLim(2);
    ylim_min(i) = sp(i).YLim(1);
    ylim_max(i) = sp(i).YLim(2);
end

xvals(1) = min(xlim_min); xvals(2) = max(xlim_max);
yvals(1) = min(ylim_min); yvals(2) = max(ylim_max);

for i = 1:numel(sp)
    set(sp(i),'XLim', xvals);
    set(sp(i),'YLim',yvals);
end
end
%% function for placing legend in extra plot
function [] = extraleg(legdat,h,nrow,ncol,s)
subplot(nrow,ncol,s+1) % initialize extra subplot and get position
legpos = get(gca,'Position');
delete( subplot(nrow,ncol,s+1)) % delete that pesky extra subplot
lg = legend(subplot(nrow,ncol,s),h,legdat);
set(lg,'Position',legpos) % stick the legend in that hole
end
%% line plotter
function [] = lnplot(txx)

yy = ylim;
% get ylims of plot for line
for sb = 1:numel(txx)
    xx = txx(sb);
    line([xx,xx],[yy(1),yy(2)],'Color',...
        [0.5,0.5,0.5,0.7]);
end
end
%% Colormap
function [cmap] = rainbow(n)
values = [
213,62,79
244,109,67
253,174,97
254,224,139
171,221,164
102,194,165
50,136,189
94,79,162]./256;

cmap = interp1(1:size(values,1), values, linspace(1,size(values,1),n), 'linear');
end
