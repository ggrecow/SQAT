function [OUTPUT_1, OUTPUT_2, OUTPUT_3] = get_EPNL(input_1, input_2, input_3, nCases)

% input parameters
method = 1;
dt = 0.5;
threshold = 10;
show = 0;

for i = 1:nCases
    
    OUTPUT_1{i} = input_1{i};
    OUTPUT_1{i}.EPNL = EPNL_FAR_Part36(input_1{i}.audio, input_1{i}.fs,... % input signal and sampling freq.
                                                            method,... % method = 0, insig is a SPL[nTime,nFreq] matrix; method = 1, insig is a sound file
                                                            dt,... % time-step in which the third-octave SPLs are averaged, in seconds.
                                                            threshold,... % threshold value used to calculate the PNLT decay from PNLTM during the calculation of the duration correction
                                                            show);
                                                        
    OUTPUT_2{i} = input_2{i};
    OUTPUT_2{i}.EPNL = EPNL_FAR_Part36(input_2{i}.audio, input_2{i}.fs,... % input signal and sampling freq.
                                                            method,... % method = 0, insig is a SPL[nTime,nFreq] matrix; method = 1, insig is a sound file
                                                            dt,... % time-step in which the third-octave SPLs are averaged, in seconds.
                                                            threshold,... % threshold value used to calculate the PNLT decay from PNLTM during the calculation of the duration correction
                                                            show);

    OUTPUT_3{i} = input_3{i};
    OUTPUT_3{i}.EPNL = EPNL_FAR_Part36(input_3{i}.audio, input_3{i}.fs,... % input signal and sampling freq.
                                                            method,... % method = 0, insig is a SPL[nTime,nFreq] matrix; method = 1, insig is a sound file
                                                            dt,... % time-step in which the third-octave SPLs are averaged, in seconds.
                                                            threshold,... % threshold value used to calculate the PNLT decay from PNLTM during the calculation of the duration correction
                                                            show);
  
end


end