PFC = []
CA1 = []
joint = []
%%
figure;hold on
% plot ((PFC./PFC(:,2))','-ob')
% plot ((CA1./CA1(:,2))','-or')
% plot ((joint./joint(:,2))','-og')


m1 = nanmean((PFC./PFC(:,2)))
m2 = nanmean((CA1./CA1(:,2)))
m3 = nanmean((joint./joint(:,2)))

e1 = nansem((PFC./PFC(:,2)))
e2 = nansem((CA1./CA1(:,2)))
e3 = nansem((joint./joint(:,2)))

errorbar([1:3],m1,e1,'-bo','LineWidth',1.5);
errorbar([1:3],m2,e2,'-ro','LineWidth',1.5);
errorbar([1:3],m3,e3,'-go','LineWidth',1.5);

ylim([0.5 1.5])
set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre-Task','Task','Post-Task'},'XTickLabelRotation',45)
xlim([0.5  3.5])

ylabel('Norm. Assembly Count')