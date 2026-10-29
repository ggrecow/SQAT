function tests = FluctuationStrength_Osses2016_reference_test
% FluctuationStrength_Osses2016 against the reference values of the four
% scripts of validation/FluctuationStrength_Osses2016 (AM tones, AM noise,
% FM tones against fmod; FM tones against the frequency deviation). As agreed
% in PR 62, the JND of 10 % is a hard limit only at the reference signal
% (unit test); here the RMSE of each case is pinned and may only improve,
% and the points inside the 10 % bars are reported. Needs the Zenodo sounds.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
tc.TestData.v = fullfile(root, 'validation');
tc.TestData.sounds = fullfile(root, 'sound_files', 'validation_SQAT_v1_0');
end

function test_rmse_against_the_references_does_not_grow(tc)
% Runs the four validation scripts of FluctuationStrength_Osses2016 (AM and FM
% tones, AM broadband noise, against the modulation frequency or the frequency
% deviation) on the Zenodo sounds and checks that the RMSE against the
% reference data stays within 2 % of its value pinned on 30.09.2026. The
% points within the 10 % JND are reported only.
tc.assumeTrue(isfolder(fullfile(tc.TestData.sounds, 'FluctuationStrength_Osses2016')), 'Zenodo sounds not found');
d = fullfile(tc.TestData.v, 'FluctuationStrength_Osses2016');
cases = {'1_AM_tones_fmod/run_validation_FS_fmod.m', 0.1252; '2_AM_BBN_fmod/run_validation_FS_AM_BBN_fmod.m', 0.2933; ...
         '3_FM_tones_fmod/run_validation_FS_FM_fmod.m', 0.7892; '4_FM_tones_freq_dev/run_validation_FS_FM_freq_dev.m', 0.1771};
for k = 1:size(cases, 1)
    S = sqat_run_script(fullfile(d, cases{k, 1}));
    tc.assertEmpty(S.run_err, S.run_err);
    if isfield(S, 'results_1khz')                 % against the frequency deviation
        got = S.results_1khz(:);
        ref = interp1(S.ref(:, 1), S.ref(:, 2), S.freq_dev(:));
    else
        got = S.results(:);
        ref = S.ref(2, :)';
    end
    r = sqrt(mean((got - ref).^2));
    name = fileparts(cases{k, 1});
    sqat_report_record('fluctuation_strength', ['rise of the RMSE over its pinned value, ' name ' (vacil)'], max(r - cases{k, 2}, 0), 0, ...
        0.02 * cases{k, 2}, 'pinned at the RMSE of 30.09.2026, 2 % for the platform (PR 62 rule)');
    inside = nnz(abs(got - ref) <= 0.1 * ref);
    sqat_report_record('fluctuation_strength', sprintf('points outside the 10 %% JND, %s (reported, %d of %d inside)', name, inside, numel(ref)), ...
        numel(ref) - inside, 0, numel(ref), 'report only: data of listening tests, not a limit of the model (PR 62)');
    tc.verifyLessThanOrEqual(r, 1.02 * cases{k, 2}, [name ': RMSE rose above its pinned value']);
end
end
