tt = 0:plot_interval:(file_info{2,File} - plot_interval);
t1 = tt(len_start:len_end);
% time_series{File} = t;
Fig = figure;
set(gcf,'outerposition',get(0,'screensize'));
for antenna = 1:n_antenna + 1
    subplot(2,2,antenna);
    hold on;
    if antenna <= n_antenna
%         plot(tt,mag_MR(:,antenna),tt,outliered_MR(:,antenna),tt,filtered_MR(:,antenna));%画出原始数据
            plot(tt,pca1_an(:,antenna));
            plot(tt,pca1_MR_an(:,antenna));
%         [~,in] = hampel(mag_MR(:,antenna),hampel_k,hampel_nsigma);
%         plot(tt(in),mag_MR(in,antenna),'sk','MarkerSize',0.2);
    else
        plot(t1, pca1);
        plot(t1, PCA_need(:,1));
        plot(t1(Idx), slopes);
        tmp = find(abs(slopes) > threshold);
        plot(t1(Idx(tmp)), slopes(tmp),'sk', 'MarkerSize',0.7);
    end
    if length(endpoint_indices{File}) > 2
        temp_odd = endpoint_indices{File}(1,1:2:size(endpoint_indices{File},2)-1);
        temp_even = endpoint_indices{File}(1,2:2:size(endpoint_indices{File},2));
        mean_sample(File) = mean(temp_even - temp_odd);
        mean_sample_length(File) = mean_sample(File) / sample_rate;
%         PCA_facial{File} = PCA_need;
        %         temp = endpoint_indices{File}(2,:) + 1;
%         temp_endpoint = endpoint_indices{File}(1,find(temp));
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
%             PCA_min = min(pca1);
%             temp_arr2 = ones(1,len) .* PCA_min;
%               stem(t(endpoint_indices{File}(1,:)), temp_arr2(endpoint_indices{File}(1,:)),...
% 'Marker','none','LineStyle','-.','Color','red');
%                     stem(t(temp_endpoint),temp_arr2(temp_endpoint),...
%             'Marker','none','LineStyle','-.','Color','green');
    end
    hold off;
    
    if antenna == n_antenna + 1
        %legend('CSI data','outlier removal data','after filter','outlier','center points');%写上对应的标签
        title_tmp = strrep([file_directory,file_info{1,File}],'\','\\');
        if antenna >= u_antenna && u_antenna ~= 0
            title(['antenna ' num2str(antenna) ...
                ';File:' strrep(title_tmp,'_','\_') ...
                ';Sample number: ' num2str(size(endpoint_indices{File},2)/2) ...
                ';Unused antenna: ' num2str(u_antenna)]);%写上标题
        else
            title(['antenna ' num2str(antenna) ...
                ';File:' strrep(title_tmp,'_','\_') ...
                ';Sample number: ' num2str(size(endpoint_indices{File},2)/2) ...
                ';Unused antenna: ' num2str(u_antenna)]);%写上标题
        end
    else
        if antenna >= u_antenna && u_antenna ~= 0
            title(['antenna ' num2str(antenna+1)]);
        else
            title(['antenna ' num2str(antenna)]);
        end
    end
    axis tight;
    xlabel('time/s');%写上对应的横纵坐标
    ylabel('amplitude');
    xlim([file_info{3,File} file_info{4,File}]);
end
set(gcf,'Color','w');%设置窗口的底纹为白色
if file_choice == 0
    if File == 1 
        delete(fullfile(file_directory, 'Figs', '*.fig'));
    end
    FigFileName = [file_info{1, File} '-' num2str(file_info{3, File}) ...
            '-' num2str(file_info{4, File})];
    savefig(Fig,fullfile(file_directory, 'Figs', num2str(File)));
    close(Fig);
end
% drawnow;
clear plot_interval t antenna in temp_* mean_* pca_max