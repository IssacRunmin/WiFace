tt = 0:plot_interval:(file_info{2,File} - plot_interval);
t1 = tt(len_start:len_end);
Fig = figure;
set(gcf,'outerposition',get(0,'screensize'));
for antenna = 1:n_antenna + 1 % 1~3 antenna; 4: PCA for all antennas
    subplot(2,2,antenna);
    hold on;
    if antenna <= n_antenna
        plot(tt,pca1_an(:,antenna));
        plot(tt,pca1_MR_an(:,antenna));
    else
        plot(t1, pca1);             % pca1 is the PCA of all antennas
        plot(t1, PCA_need(:,1));    % PCA of the Multipath removal data
        plot(t1(Idx), slopes);      % slopes are in the segmentation
        tmp = find(abs(slopes) > 10);
        plot(t1(Idx(tmp)), slopes(tmp),'sk', 'MarkerSize',0.7); % mark it
    end
    if length(endpoint_indices{File}) > 2
        temp_odd = endpoint_indices{File}...
            (1,1:2:size(endpoint_indices{File},2)-1);
        temp_even = endpoint_indices{File}...
            (1,2:2:size(endpoint_indices{File},2));
        pca_max = max(get(gca,'ylim'));
        pca_min = min(get(gca,'ylim'));
        temp_arr = ones(1,len) .* pca_max;
        
        stem(t1(temp_odd), temp_arr(temp_odd), ...
            'Marker','none','LineStyle','-.','Color','green');
        stem(t1(temp_even),temp_arr(temp_even),...
            'Marker','none','LineStyle','-.','Color','red');
        temp_arr = ones(1,len) .* pca_min;
        stem(t1(temp_odd), temp_arr(temp_odd), ...
            'Marker','none','LineStyle','-.','Color','green');
        stem(t1(temp_even),temp_arr(temp_even),...
            'Marker','none','LineStyle','-.','Color','red');
    end
    hold off;
    box on
    if antenna == n_antenna + 1
        %legend('CSI data', 'outlier removal',...
            % 'filtered','outliers','center points');
        title_tmp = strrep([file_directory,file_info{1,File}],'\','\\');
        title(['antenna ' num2str(antenna) ...
            ';File:' strrep(title_tmp,'_','\_') ...
            ';Sample number: ' num2str(size(endpoint_indices{File},2)/2)]);
    else
        title(['antenna ' num2str(antenna)]);
    end
    axis tight;
    xlabel('time/s');
    ylabel('amplitude');
    xlim([file_info{3,File} file_info{4,File}]);
end
set(gcf,'Color','w');
if file_choice == 0
    if File == 1 
        delete(fullfile(file_directory, 'Figs', '*.fig'));
    end
    FigFileName = [file_info{1, File} '-' num2str(file_info{3, File}) ...
            '-' num2str(file_info{4, File})];
    savefig(Fig,fullfile(file_directory, 'Figs', num2str(File)));
    close(Fig);
end
drawnow;
clear t antenna in temp_* mean_* pca_max