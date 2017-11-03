%Define some convenient colormaps here.
function cmap = ct_colormaps(clr)
n2 = 32;    c1 = 0.3;  c2 = 0.2;  c2b = 0;
imain = [0; linspace(c1, 1, n2)'; ones(n2, 1)];
isec =  [0; linspace(c2, 1-c2b, n2*2)'];
ilow =  [0; zeros(32, 1); linspace(0, 1-2*c2b, n2)'];

switch upper(clr)
    case {'RED','R'};               cmap = [imain, ilow, ilow];
    case {'ORANGE','O'};            cmap = [imain, isec, ilow];
    case {'YELLOW','Y'};            cmap = [imain, imain, ilow];
    case {'CHARTREUSE','CH','YG'};  cmap = [isec, imain, ilow];
    case {'GREEN','G'};             cmap = [ilow, imain, ilow];
    case {'SPRINGGREEN','SG'};      cmap = [ilow, imain, isec];
    case {'CYAN','C'};              cmap = [ilow, imain, imain];
    case {'AZURE','A'};             cmap = [ilow, isec, imain];
    case {'BLUE','B'};              cmap = [ilow, ilow, imain];
    case {'VIOLET','V'};            cmap = [isec, ilow, imain];
    case {'MAGENTA','M'};           cmap = [imain, ilow, imain];
    case {'ROSE','PINK','P'};       cmap = [imain, ilow, isec];
    case {'GRAY'};                  cmap = [isec, isec, isec];
    otherwise;                      
        try     cmap = colormap(clr);
        catch;   warning('CT_COLORMAP:badMapName', ['Invalid colormap ',...
                'or name provided.  Using default.  Type HELP CT_TRACKVIS ',...
                'or HELP COLORMAP for proper names.']);  cmap = [];
        end
end

end