function [STMtx,UnitListOut] = restrictTranges(UnitList,Trange,STmultiply)
% 
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
%
% NB - unit timestamps are assumed to be in NLX format (100us precision)
% and are convered to seconds on processing.
% Takes imported (MClust-style) .t timestamps, concatenates and pads with zeros 
%
% Dependencies: catpad.m
%
% Input arguments:
% UnitList:    Structure array of spike times, in format: data{Unit no}.t [no. spikes x 1]
% Trange:  Min and max times (s) to restict output spike trains to [min time, max time]
% STmultiply: Optionally multiply input times to bring to seconds (defaults to 1e-4 for use with mClust .t files)
%
% Output arguments:
% STMtx:   Spike time matrix with size [max no. spikes, no. units]. columns are units. shorter columns are zero-padded
% UnitListOut: Time-converted version of UnitList corrected to seconds


if nargin<2 || isempty(Trange), Trange=[0,Inf]; end; 
if nargin<3 || isempty(STmultiply), STmultiply=1e-4; end; 
UnitListOut = UnitList;
nTimeRanges = size(Trange,1);


STMtx =[]; i=1;% counter_=[];
    for iUnit=1:length(UnitList)
        temp=UnitList{iUnit}.t*STmultiply;  
        
        dataIN = []; 
        for iTRange = 1:nTimeRanges
            dataIN=[dataIN; temp(ge(temp,Trange(iTRange,1))  & lt(temp,Trange(iTRange,2)))];
        end
        dataIN=sort(dataIN);
        
        UnitListOut{iUnit}.t = dataIN./STmultiply;  % Convert spike times back to original multipliers
        
        if ~isempty(STMtx) && isnan(STMtx(1));STMtx(:,1)=[]; end % I.E. using catpad, first time round column 1 gets filled with NaNs... remove
        
        STMtx = catpad(2,STMtx,dataIN);
%         counter_=[counter_;[i length(dataIN)]];
        i=i+1;
    end
STMtx(isnan(STMtx))= 0;