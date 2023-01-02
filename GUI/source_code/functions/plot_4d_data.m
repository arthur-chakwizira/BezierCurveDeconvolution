function plot_4d_data(handles, this_data, this_slice, this_axes, file_name)
%        This function accepts the handles structure, data to plot (4D), slice
%        to plot, axes on which to plot and the file name from which the
%        data being visualised was loaded (used as title for the axes). 
%        The function generates a large number of plots from the specified
%        slice. It is called by the functions that handle loading of data
%        onto the axes. This function returns no output.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       data_size = size(this_data);
       num_plots = data_size(1)*data_size(2);
       hObject = handles.this_object;
       t = handles.fourth_dim_vector;
       this_pos = 0;
       p_bar = waitbar(0, ['Plotting curve 1 of ' num2str(num_plots)], 'Name', 'Visualing 4D data');
       plot(0, 0, 'Parent', this_axes)
       hold(this_axes, 'on')
       title(handles.axes1, [file_name ': Slice 1'], 'Interpreter', 'none')
      for x = 1:data_size(1)
            for y = 1:data_size(2) %for each pixel
                this_pos = this_pos + 1;
%                 if this_data(x,y,this_slice,:) > 100
                plot(t, squeeze(this_data(x,y,this_slice,:)),...
                'LineWidth', 0.5 , 'Parent', this_axes);
                if isgraphics(p_bar)
                    waitbar(this_pos/num_plots, p_bar, ['Plotting curve ' num2str(this_pos) ' of ' num2str(num_plots)])
                else
                    break
                end
                drawnow
%                 end
           end
      end
      if ishandle(p_bar); delete(p_bar); end
      guidata(hObject, handles)
end