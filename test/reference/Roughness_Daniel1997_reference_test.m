function tests = Roughness_Daniel1997_reference_test
% Roughness_Daniel1997 against the jury data of Daniel and Weber (1997),
% Fig. 3: AM tones of 7 carriers, run by
% validation/Roughness_Daniel1997/1_AM_modulation_freq. As agreed in PR 60,
% the JND of 17 % is a hard limit only at the reference signal (unit test);
% point by point it would reject the model itself, so the RMSE of each
% carrier is pinned and may only improve. Needs the sounds of the Zenodo
% record 7933206 in sound_files/validation_SQAT_v1_0.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
tc.TestData.v = fullfile(root, 'validation');
tc.TestData.sounds = fullfile(root, 'sound_files', 'validation_SQAT_v1_0');
end

function test_rmse_against_the_jury_does_not_grow(tc)
% Runs the validation script of Roughness_Daniel1997 (AM tones, seven
% carriers from 125 Hz to 8 kHz, against the modulation frequency) and checks
% that the RMSE of each carrier against the jury data of Daniel and Weber
% (1997) stays within 2 % of its value pinned on 30.09.2026. The points within
% the 17 % JND are reported only.
tc.assumeTrue(isfolder(fullfile(tc.TestData.sounds, 'Roughness_Daniel1997')), 'Zenodo sounds not found');
S = sqat_run_script(fullfile(tc.TestData.v, 'Roughness_Daniel1997', '1_AM_modulation_freq', 'run_validation_roughness_fmod.m'), true);
tc.assertEmpty(S.run_err, S.run_err);
il_jury(tc, S, 'roughness_daniel1997', [0.02666 0.03991 0.04944 0.02095 0.1012 0.0558 0.02857], @(r) 0.17 * r, 'the 17 % JND');
end

function il_jury(tc, S, metric, pinned, band, band_text)
% RMSE of each carrier against the jury data of Daniel and Weber (1997),
% Fig. 3, not above the pinned value (2 % allowed for the platform: the
% roughness moves by a few 1e-4 with the rounding of the FFT); the points
% inside the band are reported, not required
tags = {'125hz', '250hz', '500hz', '1khz', '2khz', '4khz', '8khz'};
x = {20:10:100, 20:10:130, 20:10:160, 20:10:160, 20:10:160, 20:10:160, 20:10:160};
r = zeros(1, 7);
inside = 0;
n = 0;
for k = 1:7
    ref = S.(['fmod_' tags{k}]);
    q = interp1(0:10:160, S.(['results_' tags{k}]), x{k}(:));
    r(k) = sqrt(mean((q - ref).^2));
    inside = inside + nnz(abs(q - ref) <= band(ref));
    n = n + numel(ref);
end
sqat_report_record(metric, 'rise of the RMSE against the jury of Daniel and Weber (1997) over its pinned value, 7 carriers (asper)', ...
    max(r - pinned, 0), zeros(1, 7), 0.02 * max(pinned), 'pinned at the RMSE of 30.09.2026, 2 % for the platform (PR 60 rule)');
sqat_report_record(metric, sprintf('points outside %s (reported, %d of %d inside)', band_text, inside, n), n - inside, 0, n, ...
    'report only: jury data, not a limit of the model (PR 60)');
tc.verifyLessThanOrEqual(r, 1.02 * pinned, 'RMSE of a carrier rose above its pinned value');
end
