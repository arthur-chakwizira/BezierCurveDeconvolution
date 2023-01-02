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

handles_update = handles;

end