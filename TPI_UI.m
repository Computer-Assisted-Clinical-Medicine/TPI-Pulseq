function TPI_UI
% ==========================================================================
% TPI_UI  Interactive front-end for the TPI Pulseq compiler.
%
% Provides parameter entry, one-click sequence compilation via
% tpi_compile_advanced(), and diagnostic visualization (full k-space
% trajectory, one-TR gradient/ADC view, and slew-vs-cone-angle sweep).
% Depends on: Pulseq (mr.*), calc_cone_distribution, design_tpi_projection,
% compute_slew_vs_cone_angle, tpi_compile_advanced.
% 
% Author:  Valentin Jost
% Date:    2026-09
% ==========================================================================

repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root, 'UI'));

fig = uifigure('Name', 'TPI Advanced Compiler', 'Position', [80 60 1500 950]);
lastOut = [];
ctl = struct();

gl = uigridlayout(fig, [2 2]);
gl.RowHeight   = {'1x', 180};
gl.ColumnWidth = {340, '1x'};

%% ---- Left: parameter panel ----
leftPanel = uipanel(gl, 'Title', 'Parameters', 'Scrollable', 'on');
leftPanel.Layout.Row = 1; leftPanel.Layout.Column = 1;

lg = uigridlayout(leftPanel, [20 2]);
lg.RowHeight = repmat({24}, 1, 20);
lg.ColumnWidth = {130, '1x'};

r = 1;
addLbl(lg, r, 'Scanner');            ctl.scanner = uidropdown(lg, 'Items', {'Skyra','Trio','TerraW60'}, 'Value', 'Skyra'); ctl.scanner.Layout.Row = r; ctl.scanner.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'Nucleus');            ctl.nuclei  = uidropdown(lg, 'Items', {'Na','H'}, 'Value', 'Na');                    ctl.nuclei.Layout.Row  = r; ctl.nuclei.Layout.Column  = 2; r = r+1;
addLbl(lg, r, 'FOV [mm]');           ctl.fov     = uieditfield(lg, 'numeric', 'Value', 240, 'Limits', [50 500]);         ctl.fov.Layout.Row     = r; ctl.fov.Layout.Column     = 2; r = r+1;
addLbl(lg, r, 'Matrix Nx');          ctl.Nx      = uieditfield(lg, 'numeric', 'Value', 60, 'Limits', [8 256], 'RoundFractionalValues', 'on'); ctl.Nx.Layout.Row = r; ctl.Nx.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'TR [ms]');            ctl.TR      = uieditfield(lg, 'numeric', 'Value', 100, 'Limits', [1 1000]);         ctl.TR.Layout.Row      = r; ctl.TR.Layout.Column      = 2; r = r+1;
addLbl(lg, r, 'Twist fraction p');   ctl.p       = uieditfield(lg, 'numeric', 'Value', 0.4, 'Limits', [0.05 1]);         ctl.p.Layout.Row       = r; ctl.p.Layout.Column       = 2; r = r+1;
addLbl(lg, r, 'T_read [ms]');        ctl.Tread   = uieditfield(lg, 'numeric', 'Value', 7.6, 'Limits', [0.1 100]);        ctl.Tread.Layout.Row   = r; ctl.Tread.Layout.Column   = 2; r = r+1;
addLbl(lg, r, 'Flip angle [deg]');   ctl.flip    = uieditfield(lg, 'numeric', 'Value', 90, 'Limits', [1 180]);           ctl.flip.Layout.Row    = r; ctl.flip.Layout.Column    = 2; r = r+1;
addLbl(lg, r, 'Ramp mode');          ctl.ramp_mode = uidropdown(lg, 'Items', {'linear','smooth'}, 'Value', 'linear');    ctl.ramp_mode.Layout.Row = r; ctl.ramp_mode.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'Ramp stretch');       ctl.ramp_stretch = uieditfield(lg, 'numeric', 'Value', 1.0, 'Limits', [1 5]);       ctl.ramp_stretch.Layout.Row = r; ctl.ramp_stretch.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'Transition mode');    ctl.transition_mode = uidropdown(lg, 'Items', {'none','taper'}, 'Value', 'none');  ctl.transition_mode.Layout.Row = r; ctl.transition_mode.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'Transition stretch'); ctl.transition_stretch = uieditfield(lg, 'numeric', 'Value', 1.0, 'Limits', [0.1 10]); ctl.transition_stretch.Layout.Row = r; ctl.transition_stretch.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'Avoid equator');      ctl.avoid_eq = uicheckbox(lg, 'Text', '', 'Value', false);                          ctl.avoid_eq.Layout.Row = r; ctl.avoid_eq.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'Symmetric cones');    ctl.sym_cones = uicheckbox(lg, 'Text', '', 'Value', true);                          ctl.sym_cones.Layout.Row = r; ctl.sym_cones.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'Undersample factor'); ctl.undersample = uieditfield(lg, 'numeric', 'Value', 1, 'Limits', [1 50]);         ctl.undersample.Layout.Row = r; ctl.undersample.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'Readout oversampling'); ctl.ro_os = uieditfield(lg, 'numeric', 'Value', 1.5, 'Limits', [1 4]);            ctl.ro_os.Layout.Row = r; ctl.ro_os.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'Echo (TEx)');         ctl.TEx = uidropdown(lg, 'Items', {'1','2'}, 'Value', '1');                         ctl.TEx.Layout.Row = r; ctl.TEx.Layout.Column = 2; r = r+1;
addLbl(lg, r, 'TE2 delay [ms]');     ctl.TE2_delay = uieditfield(lg, 'numeric', 'Value', 4.5, 'Limits', [0 50]);         ctl.TE2_delay.Layout.Row = r; ctl.TE2_delay.Layout.Column = 2; r = r+1;

r = r+1;
btnCalc = uibutton(lg, 'Text', 'Calculate / Recompile', 'FontWeight', 'bold', ...
    'BackgroundColor', [0.16 0.53 0.85], 'FontColor', 'w');
btnCalc.Layout.Row = r; btnCalc.Layout.Column = [1 2]; r = r+1;

btnSave = uibutton(lg, 'Text', 'Save Sequence (.seq)...', 'FontWeight', 'bold');
btnSave.Layout.Row = r; btnSave.Layout.Column = [1 2]; r = r+1;

%% ---- Right: Main + Diagnostics tabs ----
tg = uitabgroup(gl);
tg.Layout.Row = 1; tg.Layout.Column = 2;
tabMain = uitab(tg, 'Title', 'Main');
tabDiag = uitab(tg, 'Title', 'Diagnostics');

glMain = uigridlayout(tabMain, [1 1]);
axTraj = uiaxes(glMain); title(axTraj, 'Full 3D k-space trajectory');

glDiag = uigridlayout(tabDiag, [1 2]);
axSeq  = uiaxes(glDiag); title(axSeq, 'One TR (gradients + ADC samples)'); axSeq.Layout.Column = 1;
axSlew = uiaxes(glDiag); title(axSlew, 'Slew rate vs. cone angle');        axSlew.Layout.Column = 2;

%% ---- Bottom: status log ----
logPanel = uipanel(gl, 'Title', 'Status log');
logPanel.Layout.Row = 2; logPanel.Layout.Column = [1 2];
glLog = uigridlayout(logPanel, [1 1]);
ctl.log = uitextarea(glLog, 'Value', {'Ready.'}, 'Editable', 'off', 'FontSize', 12);

%% ---- Callbacks ----
btnCalc.ButtonPushedFcn = @(src, evt) onCalculate();
btnSave.ButtonPushedFcn = @(src, evt) onSave();

    function onCalculate()
        ctl.log.Value = {'Compiling...'};
        drawnow;

        prm.scanner             = ctl.scanner.Value;
        prm.nuclei              = ctl.nuclei.Value;
        prm.fov                 = ctl.fov.Value * 1e-3;
        prm.Nx                  = ctl.Nx.Value;
        prm.TR                  = ctl.TR.Value * 1e-3;
        prm.p                   = ctl.p.Value;
        prm.Tread               = ctl.Tread.Value * 1e-3;
        prm.flip_angle_deg      = ctl.flip.Value;
        prm.ramp_mode           = ctl.ramp_mode.Value;
        prm.ramp_stretch        = ctl.ramp_stretch.Value;
        prm.transition_mode     = ctl.transition_mode.Value;
        prm.transition_stretch  = ctl.transition_stretch.Value;
        prm.avoid_equator       = ctl.avoid_eq.Value;
        prm.symmetric_cones     = ctl.sym_cones.Value;
        prm.ro_os               = ctl.ro_os.Value;
        prm.TEx                 = str2double(ctl.TEx.Value);
        prm.TE2_delay           = ctl.TE2_delay.Value * 1e-3;
        prm.undersample_factor  = ctl.undersample.Value;

        try
            out = tpi_compile_advanced(prm);
        catch ME
            ctl.log.Value = [{'ERROR during compile:'}; {ME.message}];
            return;
        end

        ctl.log.Value = out.log(:);

        % --- Main tab: full k-space trajectory ---
        cla(axTraj); hold(axTraj, 'on');
        plot3(axTraj, out.ktraj(1,:), out.ktraj(2,:), out.ktraj(3,:), 'b-');
        plot3(axTraj, out.ktraj_adc(1,:), out.ktraj_adc(2,:), out.ktraj_adc(3,:), 'r.', 'MarkerSize', 2);
        axis(axTraj, 'equal'); grid(axTraj, 'on'); box(axTraj, 'on');
        xlabel(axTraj, 'kx'); ylabel(axTraj, 'ky'); zlabel(axTraj, 'kz');
        hold(axTraj, 'off');

        % --- Diagnostics tab: one TR view ---
        cla(axSeq); hold(axSeq, 'on');
        [wave_data, ~, ~] = out.seq.waveforms_and_times();
        colors = {'r','g','b'};
        labels = {'Gx','Gy','Gz'};
        for ch = 1:3
            if isempty(wave_data{ch}), continue; end
            t_ch = wave_data{ch}(1,:);
            a_ch = wave_data{ch}(2,:);
            mask = t_ch <= prm.TR;
            plot(axSeq, t_ch(mask)*1e3, a_ch(mask), colors{ch});
        end
        if ~isempty(out.t_adc)
            adc_mask = out.t_adc <= prm.TR;
            plot(axSeq, out.t_adc(adc_mask)*1e3, zeros(1, sum(adc_mask)), 'k.', 'MarkerSize', 4);
            labels = [labels, {'ADC samples'}];
        end
        grid(axSeq, 'on'); box(axSeq, 'on');
        xlabel(axSeq, 'time [ms]'); ylabel(axSeq, 'Gradient [Hz/m]');
        legend(axSeq, labels, 'Location', 'best');
        title(axSeq, sprintf('Gradient Waveforms', prm.TR*1e3));

        t_last_grad = 0;
        for ch = 1:3
            if ~isempty(wave_data{ch})
                t_last_grad = max(t_last_grad, max(wave_data{ch}(1, wave_data{ch}(1,:) <= prm.TR)));
            end
        end
        t_last_adc = 0;
        if ~isempty(out.t_adc)
            t_last_adc = max(out.t_adc(out.t_adc <= prm.TR));
        end
        t_active_end = max([t_last_grad, t_last_adc, 1e-6]) * 1e3 * 1.1;
        xlim(axSeq, [0 min(t_active_end, prm.TR*1e3)]);
        hold(axSeq, 'off');

        % --- Diagnostics tab: slew rate vs. cone angle ---
        cla(axSlew); hold(axSlew, 'on');
        d = out.slew_diag;
        plot(axSlew, d.theta_deg, d.Sx,          'r-',  'LineWidth', 1.0);
        plot(axSlew, d.theta_deg, d.Sy,          'g-',  'LineWidth', 1.0);
        plot(axSlew, d.theta_deg, d.Sz,          'b-',  'LineWidth', 1.0);
        plot(axSlew, d.theta_deg, d.Snorm,       'm-',  'LineWidth', 1.6);
        plot(axSlew, d.theta_deg, d.Saxis,       'k--', 'LineWidth', 1.2);
        plot(axSlew, d.theta_deg, d.Sanalytical, 'c-.', 'LineWidth', 1.4);
        plot(axSlew, d.used_theta_deg, d.used_Snorm, 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
        plot(axSlew, d.used_theta_deg, d.used_Saxis, 's', 'MarkerSize', 6, ...
            'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
        yline(axSlew, d.limit, 'k:', sprintf('Hardware limit (%.0f T/m/s)', d.limit), 'LineWidth', 1.5);
        xline(axSlew, 90, 'g:', 'kz = 0 (equator)', 'LineWidth', 1.0);
        legend(axSlew, {'S_x','S_y','S_z','S_{vector}','S_{per-axis max}', ...
            'S_{analytical} (Boada Eq.11)', 'Compiled cones (vector)', 'Compiled cones (per-axis)'}, 'Location', 'best');
        xlabel(axSlew, 'Cone polar angle \theta_0 [deg]');
        ylabel(axSlew, 'Peak slew rate [T/m/s]');
        title(axSlew, 'Slew rate vs. cone angle');
        grid(axSlew, 'on'); box(axSlew, 'on');
        hold(axSlew, 'off');

        ctl.log.Value = out.log(:);
        lastOut = out;
    end

    function onSave()
        if isempty(lastOut) || ~isfield(lastOut, 'seq')
            uialert(fig, 'Run "Calculate / Recompile" first.', 'No sequence to save');
            return;
        end
        [fn, pth] = uiputfile({'*.seq', 'Pulseq sequence file'}, 'Save TPI sequence as', 'tpi_sequence.seq');
        if isequal(fn, 0), return; end
        try
            lastOut.seq.write(fullfile(pth, fn));
            ctl.log.Value = [ctl.log.Value; {sprintf('Saved sequence to: %s', fullfile(pth, fn))}];
        catch ME
            uialert(fig, ME.message, 'Save failed');
        end
    end

end

function addLbl(parent, row, txt)
l = uilabel(parent, 'Text', txt);
l.Layout.Row = row;
l.Layout.Column = 1;
end