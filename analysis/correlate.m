function [rho, p_value] = correlate(data1, data2, mask, slice, make_plots)


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
%rho is the matrix of correlation coefficients and pval is the matrix
%of pvalues. If the two input vectors are not column vectors, corrcoef
%converts them into such. That is tmp_data1(:) and tmp_data2(:)
end