function plot_residue_functions(data, config)

tr = config.tr;
img_size = config.img_size;
slice_range = config.slice_range;
BzD = config.BzD;
oSVD = config.oSVD;
notify_when_done = config.notify_when_done;
with_dispersion = config.with_dispersion;

fitd_omega = data.fitd_omega;
r_svd = data.fitd_r_svd;
mask = data.mask;

t = 0:tr:(img_size(4)-1)*tr;


if ~BzD && ~oSVD; return; end

cprintf('*blue','PLOTTING RESIDUE FUNCTIONS \n')
cprintf('*blue', '\n')


f2 = figure(2); set(f2,'position',[200 200 550 200],'color','w');
if BzD && oSVD
    subplot(1,2,1), hold on, set(gca,'box','on','linewidth',0.1,'fontweight','bold','layer','top')
    xlabel('Time  [s]')
    ylabel('{\itR(t)}  [9]') 
    title('BzD')
    subplot(1,2,2), hold on, set(gca,'box','on','linewidth',0.1,'fontweight','bold','layer','top')
    xlabel('Time  [s]')
    ylabel('{\itR(t)}  [9]')
    title('oSVD')
    if with_dispersion
        sgtitle('With dispersion modelling')
    else
        sgtitle('No dispersion modelling')
    end
end

if oSVD && ~BzD
    figure(2); hold on, set(gca,'box','on','linewidth',0.1,'fontweight','bold','layer','top')
    xlabel('Time  [s]')
    ylabel('{\itR(t)}  [9]')
    if with_dispersion
        title('oSVD (with dispersion)')
    else
        title('oSVD (no dispersion)')
    end
end

if BzD && ~oSVD
    figure(2); hold on, set(gca,'box','on','linewidth',0.1,'fontweight','bold','layer','top')
    xlabel('Time  [s]')
    ylabel('{\itR(t)}  [9]')  
    if with_dispersion
        title('BzD (with dispersion)')
    else
        title('BzD (no dispersion)')
    end
end

tv = 0:0.1:t(end);
xrange = slice_range(1):slice_range(2);
yrange = slice_range(3):slice_range(4);
zrange = slice_range(5):slice_range(6);
% parfor_progress(img_size(1));
for x = xrange%1:img_size(1)  
    for y = yrange%1:img_size(2)
        for z = zrange %img_size_3 
            if mask(x,y,z)
                if BzD
                    r_BzD = bezier_residue_function(fitd_omega(x,y,z,:), tv);
                    figure(2)
                    if BzD && oSVD
                        subplot(1,2,1), hold on
                        plot(tv,r_BzD, '-','Color', [0.5,0.5,0.5]) 
                    else
                        plot(tv,r_BzD, '-','Color', [0.5,0.5,0.5]), hold on                        
                    end
                    axis([0 30 -0.5 1.2])
                end
                
                if oSVD                   
                    tmp_svd_r = squeeze(r_svd(x,y,z,:));    
                    tmp_svd_r = tmp_svd_r/max(tmp_svd_r); 
                    figure(2)                    
                     try
                        r_osvd = interp1(t,tmp_svd_r,tv,'pchip'); 
                        if oSVD && BzD
                            subplot(1,2,2), hold on
                            plot(tv,r_osvd,'-','Color', [0.5,0.5,0.5])
                        else
                            plot(tv,r_osvd,'-','Color', [0.5,0.5,0.5]), hold on
                        end
                     catch
                        r_osvd = tmp_svd_r;
                        if oSVD && BzD
                            subplot(1,2,2), hold on
                            plot(t,r_osvd,'-','Color', [0.5,0.5,0.5])
                        else
                            plot(t,r_osvd,'-','Color', [0.5,0.5,0.5]), hold on
                        end
                     end
                end
                
            end
        end
    end
%     parfor_progress;
end
% parfor_progress(0);

if notify_when_done
beep; beep; beep
end

cprintf('magenta','Plotting of residue functions complete. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')
end