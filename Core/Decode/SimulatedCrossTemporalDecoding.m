tb = -5:0.05:21;
times_ = num2cell(0:0.5:16);
bw = 0.05;
% sigma    = 5;    % kernel standard deviation (s) - single value
sigma    = max([1:0.05:3.2;fliplr(1:0.05:3.2)],[],1);% kernel standard deviation (s) - sweep
shoulder = 2;    % Time ranges to draw kernel over (sigmas)
plotYN = true;

templates = makeDecodingTemplates(tb,times_,bw,sigma,shoulder,plotYN);
templates_ = repmat(templates.result,1,1,10);
noise_ = abs(randn(size(templates_))*1.*rand(size(templates_)));

templates_ = templates_ +noise_;

[nCells,nT,nTrials]  = size(templates_);

LeftTrials = noise_;
RightTrials = noise_;

for i=1:nCells
    
    if rand<0.5
        LeftTrials(i,:,:) = templates_(i,:,:);
    else
        RightTrials(i,:,:) = templates_(i,:,:);
    end
    
end

figure; hold on
for i=1:length(times_)
    %     plot(tb,i*1.2+mean(squeeze(LeftTrials(i,:,:)),2),'LineWidth',1.5,'color','b')
    %     plot(tb,i*1.2+mean(squeeze(RightTrials(i,:,:)),2),'LineWidth',1.5,'color','r')
    ciplot(i*1.2+nanmean(squeeze(LeftTrials(i,:,:)),2) + nansem(squeeze(LeftTrials(i,:,:)),2),...
        i*1.2+nanmean(squeeze(LeftTrials(i,:,:)),2) - nansem(squeeze(LeftTrials(i,:,:)),2),tb,'b')
    
    ciplot(i*1.2+nanmean(squeeze(RightTrials(i,:,:)),2) + nansem(squeeze(RightTrials(i,:,:)),2),...
        i*1.2+nanmean(squeeze(RightTrials(i,:,:)),2) - nansem(squeeze(RightTrials(i,:,:)),2),tb,'r')
end
plot([0 0],[0 42],':k','LineWidth',1.5)
plot([16 16],[0 42],':k','LineWidth',1.5)
axis tight
xlabel('time (s)')
ylabel('Cell #')
%% Run statistical decoding
% FR:   Firing rate matrix [units,time x no. trials]


FR = [reshape(LeftTrials,nCells,nT*nTrials),reshape(RightTrials,nCells,nT*nTrials)]';
evt0 = [ones(1,nTrials),1+ones(1,nTrials)];
reg = 0.05;
[~,~,Ft2x,Rt2x,Ft2ciL0,Ft2ciH0,TS,dfnum,dfd]=DecodeStats(FR,evt0,reg);
TSsig = tpdf(TS,2*nTrials-2)<0.05;

figure; hold on

%     Y1=[Ft2ciL{s}(f,X1);Ft2ciH{s}(f,X1)-Ft2ciL{s}(f,X1)];
%     Y2=[Ft2ciL{s}(f,X2);Ft2ciH{s}(f,X2)-Ft2ciL{s}(f,X2)];
%     h=area(Tv(X1),Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','y'); hold on
%     h=area(Tv(X2),Y2'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','y'); hold on

h = area(tb,[Ft2ciL0;Ft2ciH0-Ft2ciL0]');set(h(1),'FaceColor','w'); set(h(2),'FaceColor','y');

plot(tb,Ft2x,'LineWidth',1.5,'color','k')
plot([0 0],[0 42],':k','LineWidth',1.5)
plot([16 16],[0 42],':k','LineWidth',1.5)
axis tight
xlabel('time (s)')
ylabel('Population decoding (F-score)')


TS(~TSsig)=0;
%     TS = smooth2a(TS,0,5);
figure; hold on
x_  = repmat(tb,size(TS,1),1);
y_  = repmat(1:size(TS,2),size(TS,1),1);
imagesc(x_(:),y_(:),TS);
set(gca,'YDir','normal')
cmap =([1 1 1;(jet)]);
colormap (cmap)
caxis([0 10])
plot([0 0],[0 size(TS,2)],':k','LineWidth',1.5)
plot([16 16],[0 size(TS,2)],':k','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Cell #')
cb = colorbar;
cb.Label.String = 'Single-cell decoding (T Score)';
axis tight
%% 
figure;
x = 0:0.1:15;
A = sum(double(TSsig))*bw;
A = cumsum(histc(A,x))./size(TSsig,2);
plot(x,A,'LineWidth',1.5,'color','k')
xlabel('Duration of significant encoding (s)')
ylabel('Fraction of cells')
set(gca,'ylim',[0,1.1])
box off
%% Run Cross temporal decoding 
% figure;plot(FR')
% [FracCorr]=DecodePredErrCrossTemporal_CV(FR,evt0,reg)
[FracCorr,FracCorr_CIshuffLow,FracCorr_CIshuffHigh]=DecodePredErrCrossTemporal_CV_GPU(FR,evt0,reg,true);
%%
decoding_ = smooth2a(FracCorr,5,5);
thresh    = smooth2a(FracCorr_CIshuffLow,5,5);
A = decoding_;
A(A<FracCorr_CIshuffLow)=NaN;
A=A>0.6;

Ticks_ =[0:4:16];



figure; hold on

imagesc(tb,tb,decoding_,'alphadata',0.3+A);
contour(tb,tb,A,'linewidth',1.5,'linecolor','k')
plot([0 0],[-5 0],'w','LineWidth',1.5);
plot([-5 0],[0 0],'w','LineWidth',1.5)
plot([max(tb)-5,max(tb)-5],[max(tb)-5 max(tb)],'w','LineWidth',1.5);
plot([max(tb)-5,max(tb)],[max(tb)-5 max(tb)-5],'w','LineWidth',1.5);
plot([-5 0],[0 0],'w')
plot([-5 max(tb)],[-5 max(tb)],'w','LineWidth',1.5)

rectangle('Position',[0 0 max(tb)-5 max(tb)-5],'EdgeColor','w','LineWidth',1.5)
rectangle('Position',[-5 -5  max(tb)+5 max(tb)+5],'EdgeColor','w','LineWidth',1.5)

scatter(0,0,50,'k')
scatter(16,16,50,'k')

set(gca,'XTick',Ticks_,'YTick',Ticks_)

%     set(gca,'XTick',Ticks_{iDelay},'YTick',Ticks_{iDelay})
% axis([-5 21 -5 21]);
colormap(jet)
%     caxis([0.4 0.6])
axis square
ylabel('Training time (s)')
xlabel('Testing time (s)')
cb = colorbar;
cb.Label.String = 'Fraction correct decoding';
axis tight



