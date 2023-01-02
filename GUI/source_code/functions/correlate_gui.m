function [rho, p_value] = correlate_gui(data1, data2, mask, slice, make_plots)
%This is a redundant function. If needed, contact
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp_data1 = data1(:,:,slice);   
tmp_data2 = data2(:,:,slice);

tmp_data1 = squeeze(tmp_data1(mask(:,:,slice))); 
tmp_data2 = squeeze(tmp_data2(mask(:,:,slice)));

if make_plots
    figure()
    plot(tmp_data1(:),tmp_data2(:),'k.')
    xlabel(num2str(inputname(1)))
    ylabel(inputname(2))
end

[rho,p_value] = corrcoef(tmp_data1,tmp_data2);
end