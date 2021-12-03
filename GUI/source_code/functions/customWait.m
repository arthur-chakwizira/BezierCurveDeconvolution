function pos = customWait(dsc_roi)
%        This function accepts a handle to a ROI object and creates a
%        listener to detect events in the ROI. Program execution waits
%        until a double click in ROI occurs, then it resumes.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hROI = dsc_roi;
% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);
% Block program execution
uiwait;
% Remove listener
delete(l);
% Return the current position
pos = hROI.Position;
end