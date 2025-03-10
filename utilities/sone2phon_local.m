function phon = sone2phon_local(sone)
% Converts loudness from sone to phon according to ISO 532-1:2017
% Input: sone - Loudness in sone
% Output: phon - Loudness level in phon
    
    % check if sone is a column vector, and correct if its not
    if size(sone,2)>1
        sone = sone';
    end

    % initialize output vector
    phon = zeros(size(sone,1),1); % column vector [nTime,1]
    
    % get idx where sone >= 1
    idx = (sone>=1);

    % calculate phon values for sone >= 1
    phon(idx) = 40 + 33.22 .* log10( sone( idx ) );

    % calculate phon values for  sone < 1
    phon(~idx) = 40 .* ( sone(~idx) + 0.0005 ).^0.35;
    
end
