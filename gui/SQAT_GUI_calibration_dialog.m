function SQAT_GUI_calibration_dialog(parent, name, c, help_text, apply)
% function SQAT_GUI_calibration_dialog(parent, name, c, help_text, apply)
%
%   Modal dialog that asks how to calibrate one sound file: full-scale level
%   (dBFS), calibrator recording or relative level (see SQAT_GUI_calibration).
%   It opens on the choice c and waits until it is closed, so that the
%   dialogs of several files come one after the other.
%
% INPUT ARGUMENTS
%   parent : the main window of the interface (the dialog opens over it)
%   name : name of the file, for the title
%   c : struct with method ('dbfs', 'calibrator', 'relative'), level (dB)
%       and file (calibrator recording), the choice the dialog opens on
%   help_text : what calibration is and the three ways, shown on top
%   apply : function handle ok = apply(method, level, calfile); OK closes
%       the dialog when it returns true, and keeps it open otherwise
%
% Cancel, or closing the window, leaves the calibration of the file as it was.
%
% Author: Sergio Aguirre and Gil Felix Greco, September 2026
%
% AI disclosure: code development in September 2026 assisted
% by Claude Opus 5 (Anthropic). All codes were verified by
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

methods = {'dbfs', 'calibrator', 'relative'};
items = {'Full-scale level (dBFS)', 'Calibrator recording', 'Relative level'};
hints = {'Level in dB SPL of a sample of value 1 (94 dB: 1 = 1 Pa).', ...
         'Level of the calibrator, e.g. 94 dB (1 Pa), 114 dB (10 Pa) or the value with an adapter.', ...
         'The rms of the whole file (all channels together) is set to this level.'};
units = {'dB SPL at full scale', 'dB SPL of the calibrator', 'dB SPL rms'};

d = uifigure('Name', ['Calibration: ' name], 'WindowStyle', 'modal', 'Tag', 'SQAT_GUI_calibration', ...
    'Position', [parent.Position(1) + 200, parent.Position(2) + 200, 560, 470], 'CreateFcn', '');
if isprop(parent, 'Theme') && ~isempty(parent.Theme)
    d.Theme = parent.Theme;
end
g = uigridlayout(d, [7 3]);
g.RowHeight = {'1x', 26, 36, 26, 26, 22, 30};
g.ColumnWidth = {110, '1x', 150};
t = uilabel(g, 'Text', help_text, 'WordWrap', 'on', 'VerticalAlignment', 'top');
t.Layout.Column = [1 3];
uilabel(g, 'Text', 'Method:', 'HorizontalAlignment', 'right');
dd = uidropdown(g, 'Items', items, 'ItemsData', methods, 'Value', c.method, 'Tag', 'cal_method', ...
    'ValueChangedFcn', @(~, ~) show());
dd.Layout.Column = [2 3];
hint = uilabel(g, 'Text', '', 'WordWrap', 'on', 'FontAngle', 'italic');
hint.Layout.Column = [1 3];
uilabel(g, 'Text', 'Level:', 'HorizontalAlignment', 'right');
lv = uieditfield(g, 'numeric', 'Value', c.level, 'Tag', 'cal_level');
unit = uilabel(g, 'Text', '');
uilabel(g, 'Text', 'Calibrator file:', 'HorizontalAlignment', 'right');
fl = uilabel(g, 'Text', c.file, 'Tag', 'cal_file', 'Tooltip', c.file);
br = uibutton(g, 'Text', 'Browse...', 'Tag', 'cal_browse', 'ButtonPushedFcn', @(~, ~) browse());
msg = uilabel(g, 'Text', '', 'FontColor', [0.76 0.12 0.17], 'Tag', 'cal_message');
msg.Layout.Column = [1 3];
uilabel(g, 'Text', '');
uibutton(g, 'Text', 'OK', 'Tag', 'cal_ok', 'ButtonPushedFcn', @(~, ~) ok());
uibutton(g, 'Text', 'Cancel', 'Tag', 'cal_cancel', 'ButtonPushedFcn', @(~, ~) delete(d));
show();
uiwait(d);

    function show()
        % the hint, the unit and the calibrator row follow the method
        k = find(strcmp(methods, dd.Value));
        hint.Text = hints{k};
        unit.Text = units{k};
        on = strcmp(dd.Value, 'calibrator');
        set([fl br], 'Enable', on);
        msg.Text = '';
    end

    function browse()
        [f, p] = uigetfile({'*.wav', 'WAV files (*.wav)'}, 'Recording of the calibrator');
        figure(d);
        if ~isequal(f, 0)
            fl.Text = fullfile(p, f);
            fl.Tooltip = fl.Text;
        end
    end

    function ok()
        if strcmp(dd.Value, 'calibrator') && isempty(fl.Text)
            msg.Text = 'Choose the recording of the calibrator.';
            return
        end
        file = '';
        if strcmp(dd.Value, 'calibrator')
            file = fl.Text;
        end
        if apply(dd.Value, lv.Value, file)
            delete(d);
        else
            msg.Text = 'The calibration could not be set: see the console output.';
        end
    end
end
