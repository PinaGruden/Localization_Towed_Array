function plot_tracks(folder,Tracks,Tracks_selected, t_serialdate, lags,parameters,data_corsscorr,data_pamguard_c,data_pamguard_w)

% Plot all and selected TDOA tracks against Pamguard detections (if 
% available) or against cross-correlograms

Nsources=size(Tracks_selected,2);

if ~isempty(folder.pamguard) % Plot against Pamguard detections
    All_data_c = data_pamguard_c;
    All_data_w =data_pamguard_w;

    figure,hold on;
    % Plot Pamguard detections:
    plot(datenum(All_data_w.time_UTC),-1.*All_data_w.tdoa,'o', ...
        'Color',[0,0,0]+0.85, 'MarkerFaceColor',[0,0,0]+0.85, 'MarkerSize', 3)
    plot(datenum(All_data_c.time_UTC),-1.*All_data_c.tdoa,'.','Color',[0,0,0]+0.85)
    % Plot All TDOA tracks
    for k=1:size(Tracks,2)
        plot(Tracks(k).time_local, Tracks(k).tdoa,'b-','LineWidth',4)
    end
    set(gca,'YDir', 'reverse');
    datetick('x','keeplimits');
    % Plot Selected TDOA tracks
    for k=1:Nsources
        plot(Tracks_selected(k).time_local, Tracks_selected(k).tdoa,'r-','LineWidth',2.5),
        text(Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1), ...
            num2str(k),'FontSize',16,'Color','k')
    end
    xlim([t_serialdate(1),t_serialdate(end)])
    xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
    h(1) = plot(NaN, NaN,'b-','LineWidth',4);
    h(2) = plot(NaN, NaN,'r-','LineWidth',2.5);
    legend(h,'All Tracked TDOAs','Selected TDOAs','Location', 'southeast');
end

if ~isempty(folder.crosscorr) % plot against cross-correlograms

    switch parameters.signal_type
        case 'both'
            Rxy_envelope_ALL_clicks=data_corsscorr{1};
            Rxy_envelope_ALL_whistles=data_corsscorr{2};

            % PLOT TRACKED TDOAS
            f=figure;
            %plot cross-correlograms (overlayed)
            ax1 = axes(f);
            im = imagesc(t_serialdate,lags,Rxy_envelope_ALL_clicks);
            datetick('x','keeplimits');
            colormap(flipud(gray(256)))
            ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
            xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
            title(['Tracked TDOAs from measurements based on ', parameters.signal_type])
            set(gca,'FontSize',14)
            caxis([0,10])
            im.AlphaData = 0.5; % change this value to change the background image transparency
            hold all;
            %plot second data
            ax2 = axes(f);
            im1 = imagesc(t_serialdate,lags,Rxy_envelope_ALL_whistles);
            datetick('x','keeplimits');
            colormap(flipud(gray(256)))
            ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
            xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
            set(gca,'FontSize',14)
            caxis([0,10])
            im1.AlphaData = 0.5; % change this value to change the foreground image transparency
            %link axes
            linkaxes([ax1,ax2])
            %%Hide the top axes
            ax2.Visible = 'off';
            ax2.XTick = [];
            ax2.YTick = [];
            set([ax1,ax2],'Position',[.17 .11 .685 .815]);
            
            %Plot All tracked TDOAs
            hold on
            for k=1:size(Tracks,2)
                plot(Tracks(k).time_local, Tracks(k).tdoa,'b-','LineWidth',4)
            end
            datetick('x','keeplimits');

            % Plot Selected TDOA tracks
            for k=1:Nsources
                plot(Tracks_selected(k).time_local, Tracks_selected(k).tdoa,'r-','LineWidth',2.5),
                text(Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1), ...
                    num2str(k),'FontSize',16,'Color','k')
            end
            xlim([t_serialdate(1),t_serialdate(end)])

            h(1) = plot(NaN, NaN,'b-','LineWidth',4);
            h(2) = plot(NaN, NaN,'r-','LineWidth',2.5);
            legend(h,'All Tracked TDOAs','Selected TDOAs','Location', 'southeast');
            set(gca,'FontSize',14)
            set(findall(gcf,'type','text'),'FontSize',14)

        otherwise

            % PLOT TRACKED TDOAS
            figure;
            %plot cross-correlogram
            imagesc(t_serialdate,lags,data_corsscorr),
            datetick('x','keeplimits');
            colormap(flipud(gray(256)))
            caxis([0,10])
            colorbar
            ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
            xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
            title(['Tracked TDOAs from measurements based on ', parameters.signal_type])

            %Plot All tracked TDOAs
            hold on
            for k=1:size(Tracks,2)
                plot(Tracks(k).time_local, Tracks(k).tdoa,'b-','LineWidth',4)
            end
            datetick('x','keeplimits');

            % Plot Selected TDOA tracks
            for k=1:Nsources
                plot(Tracks_selected(k).time_local, Tracks_selected(k).tdoa,'r-','LineWidth',2.5),
                text(Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1), ...
                    num2str(k),'FontSize',16,'Color','k')
            end
            xlim([t_serialdate(1),t_serialdate(end)])

            h(1) = plot(NaN, NaN,'b-','LineWidth',4);
            h(2) = plot(NaN, NaN,'r-','LineWidth',2.5);
            legend(h,'All Tracked TDOAs','Selected TDOAs','Location', 'southeast');
            set(gca,'FontSize',14)
            set(findall(gcf,'type','text'),'FontSize',14)

    end

end


end