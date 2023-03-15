function [ validation ] = checkACLimits(mpc, options)
%CHECKLIMITS Summary of this function goes here
%   Detailed explanation goes here
    qr = runpf(mpc, options);
    validation = 1;
    if qr.success
        if any([any(abs(qr.bus(:,8)) > qr.bus(:,12)), any(abs(qr.bus(:,8)) < qr.bus(:,13))]) && 0
            disp('The AC Power Flow voltage magnitude violation');
            validation = 0;
        elseif any(sqrt(qr.branch(:,14).^2 + qr.branch(:,15).^2) > qr.branch(:,8))
            disp('The AC Power Flow transport capacity violation');
            validation = 0;
        elseif sum(qr.gen(end/2+1:end,2)) > 1e-3
            disp('The AC Power Flow generation capacity violation');
            validation = 0;
        elseif sum(qr.gen(end/2+1:end,3)) > 1e-3
            disp('The AC Power Flow reactive generation capacity violation');
            validation = 0;
        end
    else
        validation = 0;
        disp('The AC Power Flow of the DC solution diverges ');
    end
end

