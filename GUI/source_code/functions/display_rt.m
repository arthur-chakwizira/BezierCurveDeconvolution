function display_rt(params_row, row_num, slice_num,caller, handles)
%        This function accepts a row of residue function data, row number, slice number,
%        caller (which deconvolution algorithm is running) and the handles
%        structure. It displays residue functions on the GUI axes during deconvolution.
%        This function has no outputs.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tr = handles.tr;
img_size = handles.img_size;
t = 0:tr:(img_size(4)-1)*tr;

if strcmp(caller, 'bzd')
    omega_row = squeeze(params_row(1,:,1:5));
    bzd_r_row = zeros(size(omega_row, 1), img_size(4)); 
    for r = 1:size(bzd_r_row, 1)
        bzd_r_row(r, :) = bezier_residue_function(squeeze(omega_row(r,:)), t);
        plot(t,  bzd_r_row(r, :), '-','Color', [0.5,0.5,0.5], 'LineWidth', 0.5, 'Parent', handles.rt_axes)
        title(handles.rt_axes, 'Residue functions')
        ylabel(handles.rt_axes, 'R(t)')
       xlim(handles.rt_axes, [0 50])
       ylim(handles.rt_axes, [-0.5 1.5])
       hold(handles.rt_axes, 'on')
       drawnow
    end
    if row_num == handles.slice_range(2)
        cla(handles.rt_axes, 'reset')
        set(handles.rt_axes,'XTick', [], 'YTick', [])
        set(handles.rt_axes, 'box', 'on')
    end
end

if strcmp(caller, 'svd')
    svd_r_row = squeeze(params_row(1,:,:));
    for r = 1:size(svd_r_row, 1)
        r_norm = squeeze(svd_r_row(r, :));
        plot(t,  r_norm./max(r_norm(:)), '-', 'LineWidth', 0.5, 'Parent', handles.rt_axes)
        title(handles.rt_axes, 'Residue functions')
        ylabel(handles.rt_axes, 'R(t)')
       xlim(handles.rt_axes, [0 50])
       ylim(handles.rt_axes, [-0.5 1.5])
       hold(handles.rt_axes, 'on')
       drawnow
    end    
   if row_num == handles.slice_range(2)
        cla(handles.rt_axes, 'reset')
        set(handles.rt_axes,'XTick', [], 'YTick', [])
        set(handles.rt_axes, 'box', 'on')
    end
end
