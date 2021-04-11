classdef TextureProfile
    % TextureProfile Class object containing profilometry data
    % 
    properties
        Name % Name of texture
        rHeightMap % Raw profilometry [double]
        XY_Resolution % Resolution in x & y dimensions
        
        pHeightMap % Working profilometry [double] 
    end
    
    methods
       %% Default functions
       function inspect(obj)
           
       end
        
       %% Profile Plot        
        function fig = profile_plot(obj, source, style)
        % Plots input heatmap in 2 or 3 dimensions
        % fig = profile_plot(obj, 'r'/'p', 'profile'/'3d')  
            doi = get_data(obj, source);
            
            [r_axis, c_axis] = get_axes(obj);
            set(0,'DefaultFigureWindowStyle','normal')
            
            switch style
                case 'profile'
                    fig = figure('Name', [obj.Name, ' - Profile']); 
                    set(fig, 'Position',  [100, 150, 800, 800])
                    
                    axes('Position', [0.1 0.35 0.55 0.55])
                    imagesc(r_axis, c_axis, doi);
                    display_ratio = length(c_axis) / length(r_axis);
                    pbaspect([display_ratio 1 1]); colormap jet
                    title('Height Map')
                    yticklabels({}); xticklabels({});
                    annotation('textarrow', [0.075 0.075], [0.9 0.8])
                    
                    axes('Position', [0.1 0.1 0.55 0.2])
                    plot(c_axis, mean(doi,1), 'k'); xlim([0 max(c_axis)]);
                    xlabel('mm'); ylabel('Height (um)');
                    
                    axes('Position', [0.7 0.35 0.2 0.55])
                    plot(r_axis, mean(doi,2), 'k'); xlim([0 max(r_axis)]);
                    xlab = xlabel('mm'); ylabel('Height (um)');
                    set(gca, 'View', [90 90], 'YAxisLocation', 'right','XAxisLocation', 'top')
                    set(xlab, 'Rotation', -90, 'VerticalAlignment', 'bottom')
                    
                    axes('Position', [0.7 0.1 0.2 0.2])
                    [n,b] = hist(doi(:),100);
                    plot(b,n, 'k')
                    xlabel('Height (um)'); 
                    set(gca, 'YAxisLocation', 'right')
                    ylab = ylabel('Frequency');
                    set(ylab, 'Rotation', -90, 'VerticalAlignment', 'bottom')

                case '3d'
                    fig = figure('Name', obj.Name); 
                    surf(r_axis, c_axis, doi', 'EdgeAlpha', 0);
                    daspect([1 1 1000]); set(gca, 'View', [30 20]);
                    xlim([0 max(r_axis)]); ylim([0 max(c_axis)])
                    xlabel('mm'); ylabel('mm'); zlabel('Height (um)');
                    colormap jet
                    title(obj.Name)
            end
            
        end
        
        %% Frequency Plot
        function fig = frequency_plot(obj, source)
        % Plots input heatmap in 2 or 3 dimensions
        % fig = profile_plot(obj, 'r'/'p', 'profile'/'3d')  
            set(0,'DefaultFigureWindowStyle','normal')  
            
            profile_fft = surface_fft(obj, source);
            
            fig = figure('Name', [obj.Name, ' - FFT']);
            set(fig, 'Position',  [900, 150, 400, 800])
            
            axes('Position', [0.15 0.6 0.7 0.3]); hold on
            plot(profile_fft.freq_axis, profile_fft.psd, 'k');
            title('PSD Scanning Direction'); 
            ylabel('Power'); yticklabels({})
            xlabel('Spatial Frequency (mm)'); xlim([0.25 10]); xticks([0:2:10]); 
            
            doi = get_data(obj, source);
            axes('Position', [0.15 0.1 0.7 0.4])
            temp = abs(fftshift(fft2(doi)));
            fs = 1 ./ obj.XY_Resolution;
            nyq_f = fs ./ 2;
            r_freq_ax = linspace(-nyq_f(2),nyq_f(2), size(doi,1));
            c_freq_ax = linspace(-nyq_f(1),nyq_f(1), size(doi,2));
            %[~,c_ind] = min(abs(c_freq_ax))
            
            imagesc(r_freq_ax, c_freq_ax, temp); pbaspect([1 1 1])
            title('2D FFT'); xlabel('Frequency (Hz)'); ylabel('Frequency (Hz)');
            xlim([-10 10]); ylim([-10 10])
            %caxis([0 mean(temp(:))*5])
        end
        
        %% Truncation
        function tHeightMap = h_truncate(obj, source, range)
        % Truncates input heightmap between designated range
        % tHeightMap = h_truncate(obj, 'r'/'p', [min max]) or h_truncate(obj, 'r'/'p')
            doi = get_data(obj, source);
            
            if nargin == 2
                mean_h = mean(doi(:));
                std_h = std(doi(:));
                range = [mean_h - std_h*3, mean_h + std_h*3];
                %display(['Auto truncating points below ', num2str(round(range(1),2)), ' and above ', num2str(round(range(2),2)), ' um']);
            end
            
            mask = doi < range(2) & doi > range(1);
            pre_interp = doi;
            pre_interp(mask == 0) = NaN;
            tHeightMap = fillmissing(pre_interp,'linear','EndValues','nearest');
        end
        
        
       %% Deplaning Regression Function
        function dpHeightMap = deplane(obj, source)
        % Deplanes the heightmap
        % dpHeightMap = deplane(obj, 'r'/'p')
            doi = get_data(obj, source);
            
            r_tilt = mean(doi,2)';
            c_tilt = mean(doi,1);

            p_r = polyfit([1:size(doi,1)], r_tilt,1);
            p_c = polyfit([1:size(doi,2)], c_tilt,1);
            
            corners = [0,...
                       p_c(1) * size(doi,2);... cumulative effect of column regress
                       p_r(1) * size(doi,1),... cumulative effect of row regress
                       p_r(1) * size(doi,1) + p_c(1) * size(doi,2)]; %cumulative effect of column + row regress

            plane = NaN(size(doi));
            %{
            plane(1,1) = corners(1,1); plane(1,end) = corners(1,2);
            plane(end,1) = corners(2,1); plane(end,end) = corners(2,2);
            %}
            
            plane(1,:) = linspace(corners(1,1), corners(1,2), size(plane,2));
            plane(end,:) = linspace(corners(2,1), corners(2,2), size(plane,2));
            plane = fillmissing(plane, 'linear');

            dpHeightMap = doi - plane;
        end
        
       %% Fourier Analysis
        function profile_psd = surface_fft(obj, source)
        % PSD of texture in the scanning direction
        % [profile_fft] = surface_fft(obj, 'r'/'p')
            doi = get_data(obj, source);
            
            samp_period = obj.XY_Resolution(1);
            fs = 1/samp_period; % samples per mm
            % Welch's implementation
            n_windows = 2;
            win_length = round(size(doi,1) / n_windows);
            [temp_fft, profile_psd.freq_axis] = pwelch(doi,win_length,[],[], fs);
            profile_psd.psd = mean(temp_fft,2);
            
        end
        
       %% Accessory Functions
        % Get source data
        function doi = get_data(obj, source)
            switch source
                case 'r'
                    doi = obj.rHeightMap;
                case 'p'
                    doi = obj.pHeightMap;
            end
        end
        
       % Axes values
        function [r_axis, c_axis] = get_axes(obj)
            r_axis = [0:obj.XY_Resolution(2):(size(obj.rHeightMap,1)-1)*obj.XY_Resolution(2)];
            c_axis = [0:obj.XY_Resolution(1):(size(obj.rHeightMap,2)-1)*obj.XY_Resolution(1)];            
        end
        
    end % of Methods

end % of ClassDef