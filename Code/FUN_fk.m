function fK = FUN_fk(pstar,qK,gamma)

% Author    : Abhyudaya Singh
% UIN       : 655691216
% NetID     : singh124@illinois.edu

    % First extract density, velocity, and pressure from q
    rhoK = qK(1);
    pK = qK(3);
    
    if pstar <= pK
        % Rarefaction wave
        % Calculate speed of sound
        cK = sqrt(gamma * pK / rhoK);
        % Return fK
        fK = 2 * cK / (gamma - 1) * ((pstar / pK)^((gamma - 1)/(2 * gamma)) - 1);
    else
        % Shock wave
        % Calculate A and B
        AK = 2 / (rhoK * (gamma + 1));
        BK = (gamma - 1) / (gamma + 1) * pK;
        % Return fK
        fK = (pstar - pK) * sqrt(AK / (BK + pstar));
    end
end