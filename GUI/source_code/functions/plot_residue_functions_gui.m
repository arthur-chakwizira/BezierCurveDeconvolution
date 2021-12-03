function plot_residue_functions_gui(handles)
%        This function accepts the handles structure and populates the axes
%        specified by 'rt_axes' with residue functions for a slice
%        specified by r_slice. Slider movement on the said axes will invoke
%        this function to plot residue functions for the next slise.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tr = handles.tr;
img_size = handles.img_size;
slice_range = handles.slice_range;
BzD = handles.BzD;
do_SVD = handles.do_SVD;
fitd_omega = handles.fitd_omega;
r_svd = handles.fitd_r_svd;
mask = handles.mask;
t = 0:tr:(img_size(4)-1)*tr;
wrap_text = 'Plotting residue functions...';
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')

cla(handles.rt_axes)

r_slice = handles.r_slice; %plotting residue functions for this slice only; slider will move us to next slice

tv = 0:0.1:t(end);
xrange = slice_range(1):slice_range(2);
yrange = slice_range(3):slice_range(4);
zrange = r_slice:r_slice;
for x = xrange
    for y = yrange
        for z = zrange
            if mask(x,y,z)
                if BzD
                    r_BzD = bezier_residue_function(fitd_omega(x,y,z,:), tv);
                        plot(tv, r_BzD, '-','Color', [0.5,0.5,0.5], 'LineWidth', 0.5, 'Parent', handles.rt_axes)
                        xlabel(handles.rt_axes, 'time [s]')
                        ylabel(handles.rt_axes, 'R(t)')
                        hold(handles.rt_axes, 'on')
                        drawnow
                end
                
                if do_SVD
                    tmp_svd_r = squeeze(r_svd(x,y,z,:));    
                    tmp_svd_r = tmp_svd_r/max(tmp_svd_r);
                     try

                        tmp_svd_r_fine = interp1(t,tmp_svd_r,tv,'pchip'); 
                           plot(tv, tmp_svd_r_fine, '-', 'Color', [0.5,0.5,0.5], 'LineWidth', 0.5, 'Parent', handles.rt_axes)
                            xlabel(handles.rt_axes, 'time [s]')
                            ylabel(handles.rt_axes, 'R(t)')
                           hold(handles.rt_axes, 'on')
                           drawnow
                     catch
                        tmp_svd_r_fine = tmp_svd_r;
                            plot(t, tmp_svd_r_fine, '-', 'Color', [0.5,0.5,0.5], 'LineWidth', 0.5, 'Parent', handles.rt_axes)
                            xlabel(handles.rt_axes, 'time [s]')
                            ylabel(handles.rt_axes, 'R(t)')
                           hold(handles.rt_axes, 'on')
                           drawnow
                     end
                end    
            end
        end
    end
end
wrap_text = 'Done.';
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
end