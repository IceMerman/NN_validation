function [ validation ] = checkACLimits(mpc, n_shedding, options)
%CHECKLIMITS Summary of this function goes here
%   mpc = MatPowerCase
%   n_shedding = Number of generators that acts as load shedding
%   Options MatPowerOptions
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
        elseif sum(qr.gen(end-n_shedding:end,2)) > 1e-3
            disp('The AC Power Flow generation capacity violation');
            validation = 0;
        elseif sum(qr.gen(end-n_shedding:end,3)) > 1e-3
            disp('The AC Power Flow reactive generation capacity violation');
            validation = 0;
        end
    else
        validation = 0;
        disp('The AC Power Flow of the DC solution diverges ');
    end
end

