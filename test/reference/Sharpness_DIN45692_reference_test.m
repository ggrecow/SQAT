function tests = Sharpness_DIN45692_reference_test
% Sharpness_DIN45692 against the reference values of DIN 45692:2009 for 21
% narrowband and 20 broadband noises, run by
% validation/Sharpness_DIN45692: each value within 5 % of the reference, the
% tolerance of the standard. Needs the Zenodo sounds.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
tc.TestData.v = fullfile(root, 'validation');
tc.TestData.sounds = fullfile(root, 'sound_files', 'validation_SQAT_v1_0');
end

function test_narrowband_and_broadband_noises_within_5_percent(tc)
% Runs the validation script of Sharpness_DIN45692 on the narrowband and
% broadband noises of DIN 45692:2009, Tables A.2 and A.3, and checks every
% value against its reference within 5 %.
tc.assumeTrue(isfolder(fullfile(tc.TestData.sounds, 'Sharpness_DIN45692')), 'Zenodo sounds not found');
S = sqat_run_script(fullfile(tc.TestData.v, 'Sharpness_DIN45692', 'validation_Sharpness_DIN45692_narrowband_and_broadband_signals.m'));
tc.assertEmpty(S.run_err, S.run_err);
for g = {'narrow', 'broad'}
    s = vertcat(S.(['s_' g{1}]).Sharpness)';
    r = S.(['acum_' g{1} 'band']);
    sqat_report_record('sharpness_din45692', sprintf('%sband noises against DIN 45692 (relative)', g{1}), s ./ r, ones(size(r)), 0.05, ...
        'DIN 45692:2009, 5 %');
    tc.verifyLessThanOrEqual(abs(s ./ r - 1), 0.05, [g{1} 'band noises']);
end
end
