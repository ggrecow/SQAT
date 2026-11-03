function tests = Do_SLM_reference_test
% The sound level meter of SQAT (Do_SLM, Gen_weighting_filters, Get_Leq)
% against IEC 61672-1:2013, through the four scripts of
% validation/sound_level_meter. Each script generates its own signals,
% holds its criteria and sets one ok_* flag per criterion; the tests run it
% and require every flag.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
tc.TestData.dir = fullfile(root, 'validation', 'sound_level_meter');
end

function test_frequency_weightings_a_c_and_z(tc)
% Table 3 of IEC 61672-1: the deviation of the A, C and Z filters from the
% design goal stays inside the acceptance limits of classes 1 and 2 at the
% 34 nominal frequencies (5.5.6); the transcription of Table 3 agrees with
% the expressions of Annex E within its rounding; and, as a regression
% check, the deviation stays below 0.05 dB from 10 Hz to 4 kHz.
il_require(tc, 'validation_Do_SLM_frequency_weightings.m', {'ok_table', 'ok_class1', 'ok_class2', 'ok_tight'});
end

function test_time_weighting_keeps_the_mean_square(tc)
% Exponential time weighting (IEC 61672-1, 3.7 and 5.8) is a first-order
% low-pass on the squared pressure with unit gain at DC, so the energy
% average of its output equals the mean square of the signal for any crest
% factor; the levels stay within 0.05 dB of 20 log10(rms/p0).
S = il_require(tc, 'validation_Do_SLM_time_weighting.m', {'ok'});
sqat_report_record('slm', 'time-weighted level against the rms of 7 signals of known crest factor (dB)', ...
    S.L_sqat(:), S.L_ref(:), S.tol, 'validation script: 0.05 dB from the definition (IEC 61672-1, 5.8)');
end

function test_toneburst_response(tc)
% Table 4 and Equation (7) of IEC 61672-1 (5.9): the response to single
% 4 kHz tonebursts from 0.25 ms to 1 s, time weightings F and S, A, C and Z
% weightings, stays inside the class 1 and class 2 limits, and within
% 0.2 dB of 10 log10(1 - exp(-Tb/tau)).
il_require(tc, 'validation_Do_SLM_toneburst_response.m', {'ok_class1', 'ok_class2', 'ok_tol'});
end

function test_repeated_tonebursts(tc)
% Clause 5.10 and Equation (9) of IEC 61672-1: the time-averaged level of a
% sequence of n tonebursts of duration Tb in Tm differs from the steady
% sinusoid by 10 log10(n Tb/Tm), inside the limits of Table 4 and within
% 0.2 dB of Equation (9).
il_require(tc, 'validation_Do_SLM_repeated_tonebursts.m', {'ok_class1', 'ok_class2', 'ok_tol'});
end

function S = il_require(tc, script, flags)
% runs a validation script in a workspace of its own (it clears all) and
% verifies each of its ok flags
S = sqat_run_script(fullfile(tc.TestData.dir, script));
tc.assertEmpty(S.run_err, sprintf('%s stopped: %s', script, S.run_err));
for f = flags
    tc.assertTrue(isfield(S, f{1}), sprintf('%s sets no %s', script, f{1}));
    tc.verifyTrue(S.(f{1}), sprintf('%s: %s', script, f{1}));
end
end
