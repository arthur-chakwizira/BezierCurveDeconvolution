function clickCallback(~,evt)
%This is a helper function for the function customWait.
%See customWait
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(evt.SelectionType,'double')
    uiresume
end
end