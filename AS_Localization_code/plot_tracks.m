function plot_tracks(folder,Tracks,Tracks_selected, t_serialdate, lags,parameters,tdoa_cutoff,data_corsscorr,varargin)
% plot_tracks.m is a function that plots all and selected TDOA tracks 
% against Pamguard detections (if available) or against cross-correlograms
%
% INPUTS:
% - folder: a structure specifying paths to data. At least two fields:
%           ~ folder.crosscorr - path to Cross-correlogram information.
%                               If not avaiable then folder.crosscorr = [];
%           ~ folder.pamguard - path to extracted Pamguard detections.
%                               If not avaiable then folder.pamguard = [].
% - Tracks: a structure containing all TDOA tracks. Structure has 3 fields:
%           ~ time - a vector of times (starting at 0) for a given track;
%           ~ time_local- a vector of times (serial date format) for a
%                           given track;                          
%           ~ tdoa - a vector of tdoas for a given track.
% - Tracks_selected: a structure containing selected TDOA tracks. Same fields as
%                   "Tracks".
% - t_serialdate - a vector of times (in serial date format) for the encounter. 
% - lags - a vector of all possible TDOAs (for a given sensor spacing).
% - parameters - a structure containing at least 3 fields:
%               ~ parameters.signal_type - a string specifying what signal 
%               type was used to construct cross-correlograms ('clicks'/
%               'whistles'/'both')
%               ~ parameters.d - sensor separation (in m)
%               ~ parameters.c - speed of sound (in m/s)
% - tdoa_cutoff- a scalar- TDOA cutoff value used to select tracks.
% - data_corsscorr - a 1 x N cell array containing cross-correlogram data.
%               When parameters.signal_type == 'both' -> N=2, otherwise N=1.
% - data_pamguard_c- a table containing all Pamguard detections from its
%                   click detector (if available). Table has 4 columns:
%           ~ tdoa - a single tdoa for each detected click
%           ~ bearing - a single bearing for each detected click
%           ~ time - a time (starting at 0 at the start of encounter) when
%                   click was detected
%           ~ time_UTC- a time (specified as date and time in UTC) when
%                       click was detected                      
% - data_pamguard_w- a table containing all Pamguard detections from its
%                   whistle detector (if available). Table has 4 columns
%                   and are the same as above, but for whistles.
%
%
%Pina Gruden, 2022, UH Manoa

Nsources=size(Tracks_selected,2);

b_col=[0 0.4470 0.7410];
%r_col=[0.8500 0.3250 0.0980];
y_col=[0.9290 0.6940 0.1250];
g_col=[0,0,0]+0.85;

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
            title(['Cross-correlogram with tracked and selected TDOAs from measurements based on ', parameters.signal_type])
            clim([0,10])
            im.AlphaData = 0.5; % change this value to change the background image transparency
            hold all;
            %plot second data
            ax2 = axes(f);
            im1 = imagesc(t_serialdate,lags,Rxy_envelope_ALL_whistles);
            datetick('x','keeplimits');
            colormap(flipud(gray(256)))
            ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
            xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
            clim([0,10])
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
                plot(ax2,Tracks(k).time_local, Tracks(k).tdoa,'-','Color',b_col,'LineWidth',4)
            end
            datetick('x','keeplimits');

            % Plot Selected TDOA tracks
            hold on
            for k=1:Nsources
                plot(ax2,Tracks_selected(k).time_local, Tracks_selected(k).tdoa,'-','Color',y_col,'LineWidth',2.5),
                plot(Tracks_selected(k).time_local(1), Tracks_selected(k).tdoa(1),'ko','MarkerFaceColor','y','LineWidth',1.5),
                text(ax2,Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1)-0.0005, ...
                    num2str(k),'FontSize',16,'Color','k')
            end

            % Plot the User selected TDOA cutoff:
            hold on
            plot(t_serialdate, abs(tdoa_cutoff).*ones(size(t_serialdate)),'m:','LineWidth',1.5)
            plot(t_serialdate, -1*abs(tdoa_cutoff).*ones(size(t_serialdate)),'m:','LineWidth',1.5)

            h1(1) = plot(NaN, NaN,'-','Color',b_col,'LineWidth',4);
            h1(2) = plot(NaN, NaN,'m:','LineWidth',2);
            h1(3) = plot(NaN, NaN,'-','Color',y_col,'LineWidth',2.5);
            h1(4) = plot(NaN, NaN,'ko','MarkerFaceColor','y','MarkerSize', 5,'LineWidth',2);
            legend(h1,'All Tracked TDOAs',['User selected cutoff ' ...
                '(\fontname{Courier}bearing\_cutoff \fontname{Helvetica}in ' ...
                '\fontname{Courier}specify\_parameters.m\fontname{Helvetica})'], ...
                'Selected TDOAs (based on \fontname{Courier}bearing\_cutoff\fontname{Helvetica})', ...
                'Start of selected track','Location', 'southeast');
            %             set(gca,'FontSize',14)
            %             set(findall(gcf,'type','text'),'FontSize',14)
            ss = get(0, 'Screensize'); %
            set(gcf, 'Position', [ss(1) ss(4)/2-100 ss(3) ss(4)/2]); % [left bottom width height] - move the bottom down a bit so that it's not all the way to the top.
            hold off

        otherwise

            % PLOT TRACKED TDOAS
            figure;
            %plot cross-correlogram
            imagesc(t_serialdate,lags,data_corsscorr{:}),
            datetick('x','keeplimits');
            colormap(flipud(gray(256)))
            clim([0,10])
            colorbar
            ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
            xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
            title(['Cross-correlogram with tracked and selected TDOAs from measurements based on ', parameters.signal_type])

            %Plot All tracked TDOAs
            hold on
            for k=1:size(Tracks,2)
                plot(Tracks(k).time_local, Tracks(k).tdoa,'-','Color',b_col,'LineWidth',4)
            end
            datetick('x','keeplimits');

            % Plot Selected TDOA tracks
            for k=1:Nsources
                plot(Tracks_selected(k).time_local, Tracks_selected(k).tdoa,'-','Color',y_col,'LineWidth',2.5),
                plot(Tracks_selected(k).time_local(1), Tracks_selected(k).tdoa(1),'ko','MarkerFaceColor','y','LineWidth',1.5),
                text(Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1)-0.0005, ...
                    num2str(k),'FontSize',16,'Color','k')
            end
            xlim([t_serialdate(1),t_serialdate(end)])

            % Plot the User selected TDOA cutoff:
            hold on
            plot(t_serialdate, abs(tdoa_cutoff).*ones(size(t_serialdate)),'m:','LineWidth',1.5)
            plot(t_serialdate, -1*abs(tdoa_cutoff).*ones(size(t_serialdate)),'m:','LineWidth',1.5)

            h1(1) = plot(NaN, NaN,'-','Color',b_col,'LineWidth',4);
            h1(2) = plot(NaN, NaN,'m:','LineWidth',2);
            h1(3) = plot(NaN, NaN,'-','Color',y_col,'LineWidth',2.5);
            h1(4) = plot(NaN, NaN,'ko','MarkerFaceColor','y','MarkerSize', 5,'LineWidth',2);
            legend(h1,'All Tracked TDOAs',['User selected cutoff ' ...
                '(\fontname{Courier}bearing\_cutoff \fontname{Helvetica}in ' ...
                '\fontname{Courier}specify\_parameters.m\fontname{Helvetica})'], ...
                'Selected TDOAs (based on \fontname{Courier}bearing\_cutoff\fontname{Helvetica})', ...
                'Start of selected track','Location', 'southeast');
            ss = get(0, 'Screensize'); %
            set(gcf, 'Position', [ss(1) ss(4)/2-100 ss(3) ss(4)/2]);

            hold off

    end

end


if ~isempty(folder.pamguard) % Plot against Pamguard detections
    figure,hold on;

    % Plot Pamguard detections:
    for k=1:size(varargin,2)
        plot(datenum(varargin{1,k}.time_UTC),-1.*varargin{1,k}.tdoa,'o', ...
            'Color',g_col, 'MarkerFaceColor', g_col, 'MarkerSize', 3)
    end
    
    % Plot All TDOA tracks
    for k=1:size(Tracks,2)
        plot(Tracks(k).time_local, Tracks(k).tdoa,'-','Color',b_col,'LineWidth',5)
    end
    set(gca,'YDir', 'reverse');
    datetick('x','keeplimits');
    % Plot Selected TDOA tracks
    for k=1:Nsources
        plot(Tracks_selected(k).time_local, Tracks_selected(k).tdoa,'-','Color',y_col,'LineWidth',2.5),
        plot(Tracks_selected(k).time_local(1), Tracks_selected(k).tdoa(1),'ko','MarkerFaceColor','y','LineWidth',1.5),
        text(Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1)-0.0005, ...
            num2str(k),'FontSize',16,'Color','k')
    end
    xlim([t_serialdate(1),t_serialdate(end)])
    ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
    xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
    title('Tracked and selected TDOAs along Pamguard detections.')

    % Plot the User selected TDOA cutoff:
    hold on
    plot(t_serialdate, abs(tdoa_cutoff).*ones(size(t_serialdate)),'m:','LineWidth',1.5)
    plot(t_serialdate, -1*abs(tdoa_cutoff).*ones(size(t_serialdate)),'m:','LineWidth',1.5)
  

    h1(1) = plot(NaN, NaN,'-','Color',b_col,'LineWidth',4);
    h1(2) = plot(NaN, NaN,'m:','LineWidth',2);
    h1(3) = plot(NaN, NaN,'-','Color',y_col,'LineWidth',2.5);
    h1(4) = plot(NaN, NaN,'ko','MarkerFaceColor','y','MarkerSize', 5,'LineWidth',2);
    h1(5) = plot(NaN, NaN,'o','Color',g_col, 'MarkerFaceColor',g_col, 'MarkerSize', 3);
    legend(h1,'All Tracked TDOAs',['User selected cutoff ' ...
        '(\fontname{Courier}bearing\_cutoff \fontname{Helvetica}in ' ...
        '\fontname{Courier}specify\_parameters.m\fontname{Helvetica})'], ...
        'Selected TDOAs (based on \fontname{Courier}bearing\_cutoff\fontname{Helvetica})', ...
        'Start of selected track','Pamguard detections','Location', 'southeast');
    ss = get(0, 'Screensize'); %
    set(gcf, 'Position', [ss(1) ss(4)/2-100 ss(3) ss(4)/2]);

    hold off
end

end