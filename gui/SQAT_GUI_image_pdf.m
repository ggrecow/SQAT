function SQAT_GUI_image_pdf(rgb, filename, dpi)
% function SQAT_GUI_image_pdf(rgb, filename, dpi)
%
%   Writes a PDF of one page holding the image rgb: the page has the size of
%   the image at dpi dots per inch and nothing else on it. The pixels are
%   stored without loss (Flate). The file is written here, byte by byte, so
%   that the PDF of a plot is its PNG on every machine, whatever the fonts
%   and the screen of the MATLAB that draws it.
%
% INPUT ARGUMENTS
%   rgb : [h x w x 3] uint8 image, as imread gives it
%   filename : path of the file to write, replaced when it exists
%   dpi : pixels per inch of the page (300 for the PNG of a plot at 300 dpi)
%
% Author: Sergio Aguirre and Gil Felix Greco, October 2026
%
% AI disclosure: code development in October 2026 assisted
% by Claude Fable 5.1 (Anthropic). All codes were verified by
% the authors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file is part of the SQAT toolbox and is subject
% to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt> in the SQAT
% repository root. Some files in SQAT carry a different license, always
% stated in their own header; where this file depends on them, the combined
% work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is", WITHOUT
% WARRANTY OF ANY KIND, express or implied, including but not limited to the
% warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(rgb, 3) == 1
    rgb = repmat(rgb, [1 1 3]);
end
rgb = uint8(rgb(:, :, 1:3));
[h, w, ~] = size(rgb);
raw = reshape(permute(rgb, [3 2 1]), [], 1);          % r, g, b of each pixel, the rows one after the other
if usejava('jvm')
    bo = java.io.ByteArrayOutputStream();
    ds = java.util.zip.DeflaterOutputStream(bo);
    ds.write(typecast(raw, 'int8'));
    ds.close();
    data = typecast(bo.toByteArray(), 'uint8');
    filter = '/Filter /FlateDecode ';
else
    data = raw;                                        % plain pixels: a larger file, the same page
    filter = '';
end
pw = w * 72 / dpi;                                     % the page in points
ph = h * 72 / dpi;
content = sprintf('q %.4f 0 0 %.4f 0 0 cm /Im0 Do Q', pw, ph);
fid = fopen(filename, 'w');
if fid < 0
    error('SQAT_GUI_image_pdf:open', 'Cannot write %s', filename);
end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '%%PDF-1.4\n');
offsets = zeros(1, 5);
offsets(1) = ftell(fid);
fprintf(fid, '1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n');
offsets(2) = ftell(fid);
fprintf(fid, '2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n');
offsets(3) = ftell(fid);
fprintf(fid, ['3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.4f %.4f] /Contents 4 0 R ' ...
    '/Resources << /XObject << /Im0 5 0 R >> >> >>\nendobj\n'], pw, ph);
offsets(4) = ftell(fid);
fprintf(fid, '4 0 obj\n<< /Length %d >>\nstream\n%s\nendstream\nendobj\n', numel(content), content);
offsets(5) = ftell(fid);
fprintf(fid, ['5 0 obj\n<< /Type /XObject /Subtype /Image /Width %d /Height %d /ColorSpace /DeviceRGB ' ...
    '/BitsPerComponent 8 %s/Length %d >>\nstream\n'], w, h, filter, numel(data));
fwrite(fid, data, 'uint8');
fprintf(fid, '\nendstream\nendobj\n');
xref = ftell(fid);
fprintf(fid, 'xref\n0 6\n0000000000 65535 f \n');
fprintf(fid, '%010d 00000 n \n', offsets);
fprintf(fid, 'trailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n', xref);
end
