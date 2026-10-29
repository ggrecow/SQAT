function tests = shmAuditoryFiltBank_test
% Unit tests of shmAuditoryFiltBank (utilities/ECMA418_2), the auditory
% filter bank of ECMA-418-2:2025, clause 5.1.4.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_each_band_has_0_dB_at_its_centre_frequency(tc)
% 53 bands at F(z) = (81.9289/0.1618) sinh(0.1618 z), z = 0.5 to 26.5 in
% steps of 0.5 (Formula 9). "The filter has a gain of 0 dB at the centre
% frequency F(z)", which "varies slightly for the first critical bands"
% (5.1.4.1): from z = 3.5 on the gain stays within 0.05 dB (measured
% 0.020 dB), and a tone at F(z) excites band z the most.
fs = 48000;
Fz = (81.9289/0.1618) * sinh(0.1618 * (0.5:0.5:26.5));
t = (0:fs-1)' / fs;
g = zeros(1, 53);
best = zeros(1, 53);
for z = 1:53
    y = shmAuditoryFiltBank(sin(2*pi*Fz(z)*t));
    tc.assertSize(y, [fs 53]);
    r = sqrt(mean(y(fs/2:end, :).^2)) / sqrt(0.5);       % steady state, against the input rms
    g(z) = 20 * log10(r(z));
    [~, best(z)] = max(r);
end
sqat_report_record('shm', 'filter bank gain at F(z), z = 3.5 to 26.5 (dB)', g(7:end), zeros(1, 47), 0.05, 'set by us: the gain ''varies slightly'' (ECMA-418-2, 5.1.4.1); measured 0.02 dB');
tc.verifyEqual(g(7:end), zeros(1, 47), 'AbsTol', 0.05);
tc.verifyEqual(best, 1:53);
end
