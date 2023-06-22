function [tdoa_select, selected_indx, SelectedTracks] = select_tracks(Tracks,timevec,maxtdoa)
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
for k=1:Nsources
    time_indx=ismember(timevec,Tracks(k).time);
    tdoa_measured(k,time_indx)= Tracks(k).tdoa;
end

%-------- Select which tracks you want to compute the surfaces for------
prompt = "Which track/tracks you want to localize? \n" + ...
    "Look at Figures 1 and 2 for track numbers. \n " +...
    " Note, enter a single track number or if track is fragmented, \n " + ...
    "then enter fragments as a 1 by N vector, where N is number of fragments - \n " + ...
    "i.e. enter [1,2,3] for fragments 1,2,3. \n ";

SelectedTracks =input(prompt);
selected_indx=zeros(length(SelectedTracks),Ntsteps);
for n=1:length(SelectedTracks)
ind=find(abs(Tracks(SelectedTracks(n)).tdoa)<maxtdoa);
selected_indx(n,:)=ismember(timevec,Tracks(SelectedTracks(n)).time(ind));
end
selected_indx=sum(selected_indx,1);
tdoa_select=sum(tdoa_measured(SelectedTracks,:),1,'omitnan');
tdoa_select(selected_indx==0)=NaN;

end