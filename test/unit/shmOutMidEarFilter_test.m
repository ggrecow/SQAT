function tests = shmOutMidEarFilter_test
% Unit tests of shmOutMidEarFilter (utilities/ECMA418_2), the outer and
% middle/inner ear filter of ECMA-418-2:2025, clause 5.1.3.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_response_is_the_one_of_table_1(tc)
% Table 1 of the standard gives the 8 second-order sections with 6 decimals;
% the function holds them with more digits. The magnitude response from
% 20 Hz to 20 kHz stays within 0.001 dB of the one of Table 1 (measured
% 2.4e-4 dB), in the free field (sections 1 to 8) and in the diffuse field
% (sections 3 to 8, as the text under Table 1 prescribes).
T1 = [1.015896 -1.925299  0.922118 -1.925299  0.938014
      0.958943 -1.806088  0.876439 -1.806088  0.835382
      0.961372 -1.763632  0.821788 -1.763632  0.783160
      2.225804 -1.434650 -0.498204 -1.434650  0.727599
      0.471735 -0.366092  0.244145 -0.366092 -0.284120
      0.115267  0.000000 -0.115267 -1.796003  0.805838
      0.988029 -1.912434  0.926132 -1.912434  0.914161
      1.952238  0.162320 -0.667994  0.162320  0.284244];
sos = [T1(:, 1:3) ones(8, 1) T1(:, 4:5)];
n = 2^16;
x = [1; zeros(n - 1, 1)];
f = (0:n/2)' * 48000 / n;
band = f >= 20 & f <= 20000;
for field = {'free-frontal', 1:8; 'diffuse', 3:8}'
    H = il_dB(shmOutMidEarFilter(x, field{1}), band);
    Href = il_dB(sosfilt(sos(field{2}, :), x), band);
    sqat_report_record('shm', sprintf('outer and middle ear filter, %s, against Table 1 (dB)', field{1}), H, Href, 0.001, 'set by us: Table 1 prints 6 decimals (measured 2.4e-4 dB)');
    tc.verifyEqual(H, Href, 'AbsTol', 0.001, field{1});
end
end

function H = il_dB(h, band)
% magnitude response in dB on the bins of <band>
H = 20 * log10(abs(fft(h)));
H = H(band);
end
