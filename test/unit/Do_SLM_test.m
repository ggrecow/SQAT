function tests = Do_SLM_test
% Unit tests of Do_SLM (sound_level_meter), with the acceptance limits of
% IEC 61672-1:2013 as the criteria.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_a_1kHz_tone_reads_its_level_in_A_C_and_Z(tc)
% A steady 1 kHz tone at 60 dB SPL, time weighting F, read after 1.5 s (12
% time constants). Table 3: every frequency weighting is 0 dB at 1 kHz, with
% acceptance limits of +/-0.7 dB (class 1). Clause 5.5.9: the C- and
% Z-weighted levels differ from the A-weighted level by at most 0.2 dB.
fs = 44100;
dBFS = 94;
t = (0:2*fs-1)'/fs;                                % a column, as Do_SLM expects
x = sqrt(2) * 10^((60 - dBFS)/20) * sin(2*pi*1000*t);
settled = round(1.5*fs):numel(t);
L = struct();
for w = {'A', 'C', 'Z'}
    Lw = Do_SLM(x, fs, w{1}, 'f', dBFS);
    L.(w{1}) = Lw(settled);
    what = sprintf('%s-weighted level of a 60 dB SPL tone at 1 kHz (Table 3)', w{1});
    sqat_report_record('slm', what, L.(w{1}), 60 * ones(size(L.(w{1}))), 0.7, 'IEC 61672-1:2013, Table 3, class 1 limit at 1 kHz');
    tc.verifyEqual(L.(w{1}), 60 * ones(size(L.(w{1}))), 'AbsTol', 0.7, what);
end
sqat_report_record('slm', 'C against A at 1 kHz (5.5.9)', L.C, L.A, 0.2, 'IEC 61672-1:2013, 5.5.9');
sqat_report_record('slm', 'Z against A at 1 kHz (5.5.9)', L.Z, L.A, 0.2, 'IEC 61672-1:2013, 5.5.9');
tc.verifyEqual(L.C, L.A, 'AbsTol', 0.2, 'C against A at 1 kHz');
tc.verifyEqual(L.Z, L.A, 'AbsTol', 0.2, 'Z against A at 1 kHz');
end
