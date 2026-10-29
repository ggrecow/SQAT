function tests = shmBasisLoudness_test
% Unit tests of shmBasisLoudness (utilities/ECMA418_2), the rectification,
% the rms, the nonlinearity and the threshold in quiet of ECMA-418-2:2025,
% clauses 5.1.6 to 5.1.9.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_a_block_follows_formulas_21_to_25(tc)
% A block of 100 periods of a sinusoid with rms p in band z = 8.5
% (F = 1 kHz, LTQ = 0.0086 in Table 3). Formulas 21 and 22 give a rectified
% rms of p. Formula 23, with Table 2, alpha = 1.5 and cN = 0.0211964, minus
% LTQ (Formula 25), gives the basis loudness. Footnote 9 lets cN move by
% 0.25 %; the function uses 0.0211668 * 1.00132, 0.008 % below.
fs = 48000;
Fz = (81.9289/0.1618) * sinh(0.1618 * (0.5:0.5:26.5));
v = [1 0.6602 0.0864 0.6384 0.0328 0.4068 0.2082 0.3994 0.6434];
pt = 2e-5 * 10.^((15:10:85)/20);
LTQ = 0.0086;
levels = [30 60 90];
got = zeros(size(levels));
want = zeros(size(levels));
for k = 1:numel(levels)
    p = 2e-5 * 10^(levels(k)/20);
    blk = sqrt(2) * p * sin(2*pi*1000*(0:4799)'/fs);
    [~, got(k), rmsb] = shmBasisLoudness(blk, Fz(17));
    tc.verifyEqual(rmsb, p, 'RelTol', 1e-12, 'rectified rms, Formula 22');
    want(k) = 0.0211964 * (p/2e-5) * prod((1 + (p./pt).^1.5).^(diff(v)/1.5)) - LTQ;
end
sqat_report_record('shm', 'basis loudness at 30, 60, 90 dB against Formulas 23 and 25 (rel.)', (got + LTQ)./(want + LTQ), ones(1, 3), 0.0025, 'ECMA-418-2:2025, footnote 9 (0.25 % on cN)');
tc.verifyEqual((got + LTQ)./(want + LTQ), ones(1, 3), 'AbsTol', 0.0025);
end

function test_silence_has_no_loudness(tc)
% Below the threshold in quiet the basis loudness is zero (Formula 25).
Fz = (81.9289/0.1618) * sinh(0.1618 * (0.5:0.5:26.5));
[~, N] = shmBasisLoudness(zeros(4800, 1), Fz(17));
tc.verifyEqual(N, 0);
end
