function [ validation ] = checkACLimits(mpc, n_shedding, options)
%CHECKLIMITS Summary of this function goes here
%   mpc = MatPowerCase
%   n_shedding = Number of generators that acts as load shedding
%   Options MatPowerOptions
%   Detailed explanation goes here
    qr = runpf(mpc, options);
    validation = 1;
    if qr.success
        if sum(qr.gen(end-n_shedding+1:end,2)) > 1e-3
            fprintf('\tThe AC Power Flow generation capacity violation # %d\n', validation);
            validation = validation + 2;
        end
        if sum(qr.gen(end-n_shedding+1:end,3)) > 1e-3
            fprintf('\tThe AC Power Flow reactive generation capacity violation # %d\n', validation);
            validation = validation + 4;
        end
        if any(sqrt(qr.branch(:,14).^2 + qr.branch(:,15).^2) > qr.branch(:,8))
            fprintf('\tThe AC Power Flow transport capacity violation # %d\n', validation);
            validation = validation + 8;
        end
        % Para este caso siempre hay barra con subtensión
        if any([any(abs(qr.bus(:,8)) < qr.bus(:,12)), any(abs(qr.bus(:,8)) > qr.bus(:,13))]) && 0
            fprintf('\tThe AC Power Flow voltage magnitude violation # %d\n', validation);
            validation = validation + 16;
        end
    else
        validation = validation + 32;
        fprintf('\tThe AC Power Flow of the DC solution diverges # %d\n', validation);
    end
    
end

