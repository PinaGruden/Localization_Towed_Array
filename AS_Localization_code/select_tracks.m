function [tdoa_select, selected_indx, SelectedTracks,datetime_select] = select_tracks(Tracks,timevec,maxtdoa)
% select_tracks.m allows user to select which TDOA track or which fragments
% they want to use to obtain the localization for. 
%
% INPUTS:
% - Tracks: a structure containing tdoa tracks (with fields "time", "tdoa")
% - timevec : a vector of times 
% - maxtdoa : a scalar- maximum tdoa allowed
%
% OUTPUTS:
% - tdoa_select : 1 x N vector containing tdoas of the selected source, 
%                 where N is number of time steps
% - selected_indx : 1 x N vector indicating in which time step target
%                   exists (1) and in which it does not (0), where N is 
%                   number of time steps
% - SelectedTracks : 1 x M vector of selected tracks indices with respect 
%                   to Tracks stucture, where M is number of track
%                   fragments.
%
%
% Pina Gruden, Dec 2022, UH Manoa

%---------- Re-arrange Tracks in a matrix (Nsources x Ntsteps)--------
Ntsteps=numel(timevec); %number of time steps
Nsources=size(Tracks,2);
tdoa_measured= nan(Nsources,Ntsteps);
datetime_measured=nan(Nsources,Ntsteps);
for k=1:Nsources
    t_indx_temp1=ismember(timevec,Tracks(k).time);
    t_indx_temp2=ismember(timevec,Tracks(k).time_local);

    if sum(t_indx_temp1)>0
        time_indx=t_indx_temp1;
        t_indx_flag=1;
    elseif sum(t_indx_temp2)>0
        time_indx=t_indx_temp2;
        t_indx_flag=2;
    else
        error(['Your time vector (timevec) does not match any of ' ...
            'the time stamps in your tdoa track. Check for errors and try again.'])
    end
    tdoa_measured(k,time_indx)= Tracks(k).tdoa;
    datetime_measured(k,time_indx)= Tracks(k).time_local;
end

%-------- Select which tracks you want to compute the surfaces for------
prompt = "Which track/tracks you want to localize? \n" + ...
    "Look at Figures 1 and 2 for track numbers. \n " +...
    " Note, enter a single track number or if track is fragmented, \n " + ...
    "then enter fragments as a 1 by N vector, where N is number of fragments - \n " + ...
    "i.e. enter [1,2,3] for fragments 1,2,3. \n ";

SelectedTracks =input(prompt);
%Add a ctach if user presses enter by mistake or 0
while isempty(SelectedTracks) || any(SelectedTracks==0)
prompt = "This is not a valid choice. \n" + ...
    "Selected track numbers must be bigger than 0 and not empty. Try again! \n" + ...
    "Which track/tracks you want to localize? \n" + ...
    "Look at Figures 1 and 2 for track numbers. \n " +...
    " Note, enter a single track number or if track is fragmented, \n " + ...
    "then enter fragments as a 1 by N vector, where N is number of fragments - \n " + ...
    "i.e. enter [1,2,3] for fragments 1,2,3. \n ";
SelectedTracks =input(prompt);
end
selected_indx=zeros(length(SelectedTracks),Ntsteps);
for n=1:length(SelectedTracks)
    ind=find(abs(Tracks(SelectedTracks(n)).tdoa)<maxtdoa);
    if t_indx_flag==1
        selected_indx(n,:)=ismember(timevec,Tracks(SelectedTracks(n)).time(ind));
    elseif t_indx_flag==2
        selected_indx(n,:)=ismember(timevec,Tracks(SelectedTracks(n)).time_local(ind));
    end
end
selected_indx=sum(selected_indx,1);
tdoa_select=sum(tdoa_measured(SelectedTracks,:),1,'omitnan');
tdoa_select(selected_indx==0)=NaN;
datetime_select=sum(datetime_measured(SelectedTracks,:),1,'omitnan');
datetime_select(selected_indx==0)=NaN;

end