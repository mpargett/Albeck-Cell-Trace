%CT_GETPHASE
%

function out = ct_getphase(tf_est)
%   Ensure cell encapsulation
if isnumeric(tf_est); tf_est = {tf_est}; end

out = cell(size(tf_est));   %Initialize
for sp = 1:numel(tf_est)
    %Evaluate phase in degrees
    phs = 180./pi().*( atan2( imag(tf_est{sp}), real(tf_est{sp}) ) );
    %Stack discontinuities to evaluate phase past 360 degrees
    dphs = phs(2:end) - phs(1:end-1); discon = find(abs(dphs) > 180);
    for ss = 1:numel(discon)
        phs(1+discon(ss):end) = phs(1+discon(ss):end) ...
            - sign(dphs(discon(ss))).*360;
    end
    out{sp} = phs;
end
end