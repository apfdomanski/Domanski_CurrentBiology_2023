%% FINAL - Plot all activations on continuous timesale with units overlaid - one example assembly showing non-member units 
SF = 1.2;
SF_Ass = 5;
bounds_ = 10;
Offset = 4;
plotNonMembers = false;
tExtreme = [5240 5490];
% tExtreme = [5000 8000];
shading_ = [0.3 0.3];
TrialSkip = 2;
tlims_ = [...
        t.Short.CueLight_LeftCorrect,  t.Short.SamplePress_LeftCorrect,    t.Short.ChoicePress_LeftCorrect;...
        t.Short.CueLight_RightCorrect, t.Short.SamplePress_RightCorrect,   t.Short.ChoicePress_RightCorrect;...
        t.Short.CueLight_LeftError',   t.Short.SamplePress_LeftError',     t.Short.ChoicePress_LeftError';...
        t.Short.CueLight_RightError',  t.Short.SamplePress_RightError',    t.Short.ChoicePress_RightError';...
        t.Medium.CueLight_LeftCorrect, t.Medium.SamplePress_LeftCorrect,   t.Medium.ChoicePress_LeftCorrect;...
        t.Medium.CueLight_RightCorrect,t.Medium.SamplePress_RightCorrect,  t.Medium.ChoicePress_RightCorrect;...
        t.Medium.CueLight_LeftError',  t.Medium.SamplePress_LeftError',    t.Medium.ChoicePress_LeftError';...
        t.Medium.CueLight_RightError', t.Medium.SamplePress_RightError',   t.Medium.ChoicePress_RightError';...
        t.Long.CueLight_LeftCorrect,   t.Long.SamplePress_LeftCorrect,     t.Long.ChoicePress_LeftCorrect;...
        t.Long.CueLight_RightCorrect,  t.Long.SamplePress_RightCorrect,    t.Long.ChoicePress_RightCorrect;...
        t.Long.CueLight_LeftError',    t.Long.SamplePress_LeftError',      t.Long.ChoicePress_LeftError';...
        t.Long.CueLight_RightError',   t.Long.SamplePress_RightError',     t.Long.ChoicePress_RightError']*1e-6;
    
[tlims_ ,I] = sortrows(tlims_,1);
LR_         = t.LR(I');

LR_(min(tlims_,[],2)<tExtreme(1) | max(tlims_,[],2)>tExtreme(2)) = [];
tlims_(tlims_<tExtreme(1) | tlims_>tExtreme(2)) = NaN;
tlims_(sum(isnan(tlims_),2)>0,:)=[];



idx_ = closest(Tmtx,tExtreme);
Tmtx_ = Tmtx(idx_(1):idx_(2));  


AssCol = {[0.081 0.09 0.631],[0.761 0.505 0.046],[0.2,0.67,0.4]};
% [0.16,0.36,0.19]
% [0.2,0.67,0.4]    

figure('color','w'); hold on

 %%%%%%%%%%%% Plot the units

        
 iFR_ = iFR{3};%(:,units_Real);
 iFR_ = zscore(iFR_);
 
 iFR_ = bsxfun(@minus,iFR_,iFR_(1,:));
 iFR_ = iFR_(idx_(1):idx_(2),:);
 
%  units_ = Ass.units{iArea}{iAss};
 %
 
 units_Real = Ass.usel_out{3};
 for iUnit =1:length(units_Real)
     temp_ = mat2gray(iFR_(:,units_Real(iUnit)));
     temp_([1,end]) = 0;%min(temp_([1,end]));
     patch(Tmtx_,Offset+temp_+iUnit*SF,[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6],'FaceAlpha',0.6)
 end

%%
AssCollapsed = cell2mat(Ass.FSCsel);
for iArea = 3%1:3;
 
for iAss=1:size(AssCollapsed,2)

    %%%%%%%%%%%% plot the assemblies
    for iTrial =1:TrialSkip:size(tlims_,1)
        try
            idx=[];
            try
                idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
            catch
                idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
                idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
            end
            idx   = idx(1):idx(2);
            temp_ = AssCollapsed(:,iAss);
            temp_ = mat2gray(temp_);
            temp_ = temp_(idx);
            temp_([1,end]) =0;%min(temp_([1,end]));
            
            patch(Ass.Tmtx(idx),SF_Ass*(temp_)+iAss*SF_Ass - size(AssCollapsed,2)*SF_Ass,AssCol{iArea},'EdgeColor',AssCol{iArea},'FaceAlpha',0.6)
            
                  
        end
    end
    
    
%    
% 
% %%%%%%%%%%%% plot the lever presses
% for iTrial =1:TrialSkip:size(tlims_,1) 
%     try
%     idx=[];
%     try
%         idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
%     catch
%         idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
%         idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
%     end
%     idx = idx(1):idx(2);      
%     
%  
%     plot([tlims_(iTrial,2);tlims_(iTrial,2)],[-1 -0.2]*SF,'color',[0 1 0 0.6],'LineWidth',2)
%     plot([tlims_(iTrial,3);tlims_(iTrial,3)],[-1 -0.2]*SF,'color',[1 0 0 0.6],'LineWidth',2)
%     switch LR_(iTrial)
%         case 1
%             text(mean(tlims_(iTrial,2:3)),-2,'L','HorizontalAlignment','center')
%         case 2
%             text(mean(tlims_(iTrial,2:3)),-2,'R','HorizontalAlignment','center')
%     end
%     
%     end
% end
% % max_ = SF*size(Ass.FSCsel{iArea},2);

% set(gca,'xtick',[min(tlims_)-bounds_:60:max(tlims_)+bounds_])
% grid on
% axis([tExtreme -2 Inf])
% axis off
end
end

