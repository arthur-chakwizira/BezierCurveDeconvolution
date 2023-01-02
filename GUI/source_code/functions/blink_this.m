function blink_this(text, colour)
%       This function blinks the text displayed on the GUI main window.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for blink = 1:3
    set(text, 'ForegroundColor', 'w')
    pause(0.13)
    set(text, 'ForegroundColor', colour)
    pause(0.13)
end
end