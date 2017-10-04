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
%
% -------------------------------------------------------------------------
% VERSION 1.0
function [] = plot_by(plotby,data,channel,pmd)
%% Check input
plotby = lower(plotby);

if any(strcmp(plotby,{'treatment','treatments','tx','txs'}))
    plotby = 'treatment';
elseif any(strcmp(plotby,{'cell','celltype','celltypes'}))
    plotby = 'celltype';
else error('Accepted plotby types: treatment or celltype. plotby must be a string')
end

%% get cell type names
c1 = pmd.Cell(:,:,1); 
c2 = c1(cellfun(@isstr,c1)); celltypes = unique(c2); % find unique names of celltypes
celltypes = celltypes(~cellfun(@isempty,celltypes)); % discard empty name fields
% get rid of @ density if present for fieldname
cellfn = cellfun(@(x)x{1},regexp(celltypes,'@','split'),'un',0); 
% remove any spaces in name
cellfn = arrayfun(@(x)strrep(cellfn{x},' ',''),1:size(cellfn),'un',0);
% make list of xys for each celltype
for s = 1:numel(cellfn)
    %   Assemble list of xys with the current cname
    idx.Cell.(cellfn{s}) = [pmd.xy{strcmp(celltypes{s},c1)}];  
end

%% get treatments and find unique combos
catTxnames = {'pTx','Tx1','Tx2','Tx3'}; % possible treatment fields
[~,I]= intersect(catTxnames, fieldnames(pmd)); % tx fields in experiment
I = sort(I); catTxnames = catTxnames(I);
Txcat = cell(size(pmd.Cell,1),size(pmd.Cell,2)); % initialize cell array 

for s = 1:numel(catTxnames) % cat each treatment in a loop
   Txcat = cat(3,Txcat, cellfun(@(x)cat(2,x{:}),...
        cat(3,num2cell(pmd.(catTxnames{s})(:,:,1:4),3)),'Un',0));
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

% find the unique treatment combos
txs = unique(Txcat2(cellfun(@isstr,Txcat2)));
% remove extra spaces in name caused by NaNs to make pretty names
treatments = arrayfun(@(x)regexprep(txs{x},char(0),''),1:size(txs),'un',0)';
% remove double +s 
treatments = arrayfun(@(x)regexprep(treatments{x},'++','+'),1:size(treatments),'un',0)';

% vector of xys for indexing during plotter
xymat = cell2mat(pmd.xy(stridx)); 
%% Plotting
switch plotby
% Plot by Cell type
    case 'celltype'
        [nrow,ncol,pos] = subplot_sizer(numel(cellfn)); % optimal subplot dimensions
        cmap = hsv(numel(txs));
        figure,
        for s = 1:numel(cellfn)
            subplot(nrow,ncol,s), hold on
            sc = 0; legvec = [];
            for sa = 1:numel(txs) % plot line for each unique treatment
                xy =xymat(strcmp(txs{sa},Txcat2(stridx))); % find xy
                xy = intersect(xy, idx.Cell.(cellfn{s})); % find xy for treatment/celltype
                for sb = 1:length(xy) % plot line for each xy
                    sc = sc+1; % index for saving individual line objects for input to legend.
                    if sb == 1
                        legvec = [legvec,1]; % identify unique legend entries.
                    else
                        legvec = [legvec,0];
                    end
                   h{sc} = plot(nanmean(data{xy(sb)}.data.(channel)),'Color',cmap(sa,:));
                end
            end
%             xlim([0,size(data{xy(sb)}.data.(channel),2)])
            axis tight
            title(cellfn{s})
        end
        legvec = logical(legvec);
        subplot_standardizer(gcf) % make axis limits the same for all plots
        set(gcf,'Position',pos) % figure bigger
        extraleg(treatments,[h{legvec}],nrow,ncol,s) % make legend in extra subplot
        
% Plot by treatments
    case 'treatment'
        [nrow,ncol,pos] = subplot_sizer(numel(txs)); % optimal subplot dimensions
        cmap = lines(numel(cellfn)+1);
        figure,
        for s = 1:numel(txs)
            subplot(nrow,ncol,s), hold on
            sc = 0; legvec = [];
            for sa = 1:numel(cellfn) % plot line for each unique treatment
                xy = xymat(strcmp(txs{s},Txcat2(stridx))); % find xy
                xy = intersect(xy, idx.Cell.(cellfn{sa})); % find xy for treatment/celltype
                for sb = 1:length(xy) % plot line for each xy
                     sc = sc+1; % index for saving individual line objects for input to legend.
                    if sb == 1
                        legvec = [legvec,1]; % identify unique legend entries.
                    else
                        legvec = [legvec,0];
                    end
                   h{sc} = plot(nanmean(data{xy(sb)}.data.(channel)),'Color',cmap(sa,:));
                end
            end
%             xlim([0,size(data{xy(sb)}.data.(channel),2)])
            axis tight
            title(treatments{s})
        end
        legvec = logical(legvec);
        subplot_standardizer(gcf) % make axis limits the same for all plots
        set(gcf,'Position',pos) % figure bigger
        extraleg(cellfn,[h{legvec}],nrow,ncol,s) % make legend in extra subplot

end

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
