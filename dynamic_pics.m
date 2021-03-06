function dynamic_pics(samples, name_file, MU_ALL, SIG_ALL, THETA_ALL, SW_ALL, LL_ALL)
    draw_num = length(LL_ALL);
    mapp = colormap('jet');
    [K,~] = size(MU_ALL{1}{1});
    vis_num = length(MU_ALL);
    v = VideoWriter('record.avi');
    v.FrameRate = 5;
    open(v)
    f = figure(2);
    set(gcf,'position',[0,0,1920,468]);    
    maxite = length(MU_ALL{1});
    for i = 1:20:maxite
        subplot(1,3,1)
        plot(samples(:,1), samples(:,2),'.','color','b'); hold on    
        for j = 1:vis_num
            draw_direct_adv(samples, THETA_ALL{j}{i}, mapp(floor(j*64/draw_num),:))
            for k = 1:K
                error_ellipse_adv(SIG_ALL{j}{i}(1:2,1:2,k), MU_ALL{j}{i}(k,1:2), 0.95, 'style', mapp(floor(j*64/draw_num),:)); hold on
            end                
        end
        hold off
        xlim([-max(max(samples(:,1:2))) max(max(samples(:,1:2)))]);ylim([-max(max(samples(:,1:2))) max(max(samples(:,1:2)))])
        set(gca, 'fontsize', 20, 'fontname', 'Times New Roman', 'fontweight', 'Bold');
        drawnow

        subplot(1,3,2);
        for j = 1:vis_num
            semilogy(LL_ALL{j}(1:i), 'color', mapp(floor(j*64/draw_num),:), 'LineWidth', 5); hold on
        end
        hold off
        title('Likelihood v.s. Iterations')
        legend(name_file);
        set(gca, 'fontsize', 20, 'fontname', 'Times New Roman', 'fontweight', 'Bold');
        drawnow

        subplot(1,3,3);
        for j = 1:vis_num
            semilogy(SW_ALL{j}(1:i), 'color', mapp(floor(j*64/draw_num),:), 'LineWidth', 5); hold on
        end
        hold off
        title('Wasserstein v.s. Iterations')
        legend(name_file);
        set(gca, 'fontsize', 20, 'fontname', 'Times New Roman', 'fontweight', 'Bold');  
        drawnow
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    close(v)