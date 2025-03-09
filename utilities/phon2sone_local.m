function sone = phon2sone_local(phon)
% Converts loudness from phon to sone according to ISO 532-1:2017
% Output: phon - Loudness level in phon
% Input: sone - Loudness in sone

    % check if phon is a column vector, and correct if its not
    if size(phon,2)>1
        phon = phon';
    end

    % initialize output vector
    sone = zeros(size(phon,1),1); % column vector [nTime,1]
    
    % get idx where phon >= 40
    idx = (phon>=40);

    % calculate sone values for phon >= 40
    sone(idx) =  2.^(0.1*( phon(idx) - 40 ) );  

    % calculate phon values for phon < 40
    sone(~idx) = ( phon( ~idx ) ./ 40 ).^(1/0.35);

end
