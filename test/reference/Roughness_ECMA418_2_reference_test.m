function tests = Roughness_ECMA418_2_reference_test
% Roughness_ECMA418_2 against the jury data of Daniel and Weber (1997),
% Fig. 3, run by validation/Roughness_ECMA418_2/1_AM_modulation_freq. The
% figure there draws a band of 0.1 asper; as for Roughness_Daniel1997 (PR 60),
% the RMSE of each carrier is pinned and may only improve, and the points
% inside the band are reported. Needs the Zenodo sounds.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
tc.TestData.v = fullfile(root, 'validation');
tc.TestData.sounds = fullfile(root, 'sound_files', 'validation_SQAT_v1_0');
end

function test_rmse_against_the_jury_does_not_grow(tc)
% Runs the validation script of Roughness_ECMA418_2 (AM tones, seven carriers
% from 125 Hz to 8 kHz, against the modulation frequency) and checks that the
% RMSE of each carrier against the jury data of Daniel and Weber (1997) stays
% within 2 % of its value pinned on 30.09.2026. The points within 0.1 asper
% are reported only.
tc.assumeTrue(isfolder(fullfile(tc.TestData.sounds, 'Roughness_Daniel1997')), 'Zenodo sounds not found');
S = sqat_run_script(fullfile(tc.TestData.v, 'Roughness_ECMA418_2', '1_AM_modulation_freq', 'run_validation_roughness_fmod.m'), true);
tc.assertEmpty(S.run_err, S.run_err);
il_jury(tc, S, 'roughness_ecma418_2', [0.01228 0.01283 0.03765 0.06274 0.04975 0.04446 0.04634], @(r) 0.1 + 0*r, 'the band of 0.1 asper');
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
