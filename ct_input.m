%CT_INPUT
%   Input parsing for Cell Trace functions
%   
%Usage:
%   P = ct_input(A,P)
%       returns the modified parameter structure P with entries added (or
%       altered) by the Name/Value pairs in cell array A.

function p = ct_input(a,p)
%Parse inputs
nin = length(a);     %Check for even number of add'l inputs
if nin == 1 && isstruct(a)
    aa = a;  a(1:2:nin) = fieldnames(aa);  a(2:2:nin) = struct2cell(aa);
elseif rem(nin,2) == 1
    warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  
end

%Assign input pairs to parameter structure
for s = 1:2:nin;   p.(lower(a{s})) = a{s+1};   end

end