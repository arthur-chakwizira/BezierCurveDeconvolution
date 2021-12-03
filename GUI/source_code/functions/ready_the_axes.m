function handles_update = ready_the_axes(handles)
%        When deployed, this function prepares the axes for the display of
%        results, depending on what results the user wants to display. The
%        input to this function is the handles structure, and the output is
%        an updated version of the same.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%turn on all relevant image-related objects
set(handles.slider1, 'Visible', 'on')%slice change on axes1
set(handles.slider2, 'Visible', 'on')%slice change on axes2
set(handles.slider6, 'Visible', 'on') %upper window setting on axes1
set(handles.slider10, 'Visible', 'on')%lower window setting on axes1
set(handles.slider7, 'Visible', 'on')%upper window setting on axes2
set(handles.slider11, 'Visible', 'on')%lower window setting on axes2
set(handles.slider8, 'Visible', 'off')%horizontal scrolling on axes1
set(handles.slider9, 'Visible', 'off')%horizontal scrolling on axes2
set(handles.pushbutton15, 'Visible', 'on')%colormap on axes1
set(handles.pushbutton16, 'Visible', 'on')%colormap on axes2
set(handles.togglebutton2, 'Visible', 'on')%colorbar on axes1
set(handles.togglebutton3, 'Visible', 'on')%colorbar on axes2
set(handles.pushbutton17, 'Visible', 'on')%rotate on axes1
set(handles.pushbutton18, 'Visible', 'on')%rotate on axes2
%__________________________________________________________________________
%intial colour limits (window settings); can be changed using sliders
climcbf   = [0 1000];
climmtt = [0 15];
climdelay = [0 10];
clim_oef = [0 0.05];
clim_cmro2 = [0 0.5];
clim_r10 = [0 30];
clim_r50 = [0 10];
%__________________________________________________________________________
%initialise images with zeros in all the following cases
if handles.show_cbf
    handles.cbf_image = imshow(zeros(handles.img_size(1),handles.img_size(2)), climcbf, 'colormap', gray, 'Parent', handles.cbf_axes); 
    title(handles.cbf_axes, 'CBF [ml/100g/min]')
end

if handles.show_mtt
    handles.mtt_image = imshow(zeros(handles.img_size(1),handles.img_size(2)), climmtt, 'colormap', gray, 'Parent', handles.mtt_axes); 
    title(handles.mtt_axes, 'MTT [s]')
end

if handles.show_delay
    handles.delay_image = imshow(zeros(handles.img_size(1),handles.img_size(2)), climdelay, 'colormap', gray, 'Parent', handles.delay_axes); 
    title(handles.delay_axes, 'Delay [s]')
end

if handles.show_oef
  
    handles.oef_image = imshow(zeros(handles.img_size(1),handles.img_size(2)), clim_oef, 'colormap',gray, 'Parent', handles.oef_axes); 
    title(handles.oef_axes, 'OEF')
end

if handles.show_cmro2
    handles.cmro2_image = imshow(zeros(handles.img_size(1),handles.img_size(2)), clim_cmro2, 'colormap', gray, 'Parent', handles.cmro2_axes); 
    title(handles.cmro2_axes, 'CMRO2 [ml/100g/min]')
end

if handles.show_r10
    handles.r10_image = imshow(zeros(handles.img_size(1),handles.img_size(2)), clim_r10, 'colormap', gray, 'Parent', handles.r10_axes); 
    title(handles.r10_axes, 'R10 [s]')
end

if handles.show_r50
    handles.r50_image = imshow(zeros(handles.img_size(1),handles.img_size(2)), clim_r50, 'colormap', gray, 'Parent', handles.r50_axes); 
    title(handles.r50_axes, 'R50 [s]')
end
%__________________________________________________________________________

if handles.plot_residue_funcs %if user is plotting residue functions
    if handles.rt_axes == handles.axes1 %turn off all image-related objects on chosen axes
        set(handles.slider1, 'Visible', 'on')%slice change on axes1
        set(handles.slider6, 'Visible', 'off') %upper window setting on axes1
        set(handles.slider10, 'Visible', 'off')%lower window setting on axes1
        set(handles.slider8, 'Visible', 'off')%horizontal scrolling on axes1
        set(handles.pushbutton15, 'Visible', 'off')%colormap on axes1
        set(handles.pushbutton17, 'Visible', 'off')%rotate on axes1
        set(handles.togglebutton2, 'Visible', 'off')%colorbar on axes1
    end
    if handles.rt_axes == handles.axes2 %
        set(handles.slider2, 'Visible', 'on')%slice change on axes2
        set(handles.slider7, 'Visible', 'off') %upper window setting on axes2
        set(handles.slider11, 'Visible', 'off')%lower window setting on axes2
        set(handles.slider9, 'Visible', 'off')%horizontal scrolling on axes2
        set(handles.pushbutton16, 'Visible', 'off')%colormap on axes2
        set(handles.pushbutton18, 'Visible', 'off')%rotate on axes2
        set(handles.togglebutton3, 'Visible', 'off')%colorbar on axes2        
    end
    plot(0, 0, 'Parent', handles.rt_axes) %initialise plots with [0, 0]
    title(handles.rt_axes, 'Residue functions')
    xlabel(handles.rt_axes, 'time [s]')
    ylabel(handles.rt_axes, 'R(t)')
    xlim(handles.rt_axes, [0 50])
    ylim(handles.rt_axes, [-1.5 1.5])
end
handles_update = handles;

end