% clear workspace
clear all, close all, clc
cd('/Users/christopherpiecuch/Documents/WHOI_Daily_Notes_Laptop/2023/2023 06/SummerSchoolGIA/code/')

% load data
load 20230628_data4.mat % sea-level reconstructions from Engelhart and Horton (2012) in Southern Massachusetts
load clr0.mat % colormap
clr0(8,:)=mean(clr0,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1 data and location
fig=figure('color','white')
fig.Position(3) = fig.Position(3)*4/3;
fig.Position(4) = fig.Position(4)*2/3;

subplot(1,2,1)
 originLat = dm2degrees([30 0]);
 originLon = dm2degrees([-80 0]);
 axesm ('ortho', 'Frame', 'on', 'Grid', 'on','Origin',[originLat originLon]);
 geoshow('landareas.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])
 axis off
 plotm(42,-70,'p','markerfacecolor',clr0(4,:)*.5+[.5 .5 .5],'markeredgecolor',clr0(4,:),'markersize',15)
 %legend(['Study location'],'location','southwest')

subplot(1,2,2)
 box on, hold on, grid on
 plot([100 200],[100 200],'color',clr0(7,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(6,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(5,:),'linewidth',2)

 ii=[]; ii=find(ind==0);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 c=fill([x1 x2 x2 x1 x1]*1e-3,[y1 y1 y2 y2 y1],'k','linewidth',2);
 set(c,'EdgeColor',clr0(7,:),'facecolor',clr0(7,:))
 end

 % marine
 ii=[]; ii=find(ind==-1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y2 y2],'color',clr0(6,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(6,:),'linewidth',2)
 end

 % terrestrial
 ii=[]; ii=find(ind==1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y1 y1],'color',clr0(5,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(5,:),'linewidth',2)
 end
 axis([0 12 -40 0])
 alpha(0.1)

 legend([{'Index point'};{'Marine limit'};{'Freshwater limit'}],'location','south','orientation','vertical'), %legend boxoff
 xlabel('Time (kyr BP)','fontsize',12,'fontweight','normal')
 ylabel('RSL (m)','fontsize',12,'fontweight','normal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 2 ordinary least squares 
fig=figure('color','white')
fig.Position(3) = fig.Position(3)*4/3;
fig.Position(4) = fig.Position(4)*2/3;

subplot(1,2,2)

 box on, hold on, grid on
 plot([100 200],[100 200],'color',clr0(7,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(6,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(5,:),'linewidth',2)
 plot([100 200],[100 200],'color',[.5 .5 .5],'linewidth',2)

 ii=[]; ii=find(ind==0);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 c=fill([x1 x2 x2 x1 x1]*1e-3,[y1 y1 y2 y2 y1],'k','linewidth',2);
 set(c,'EdgeColor',clr0(7,:),'facecolor',clr0(7,:))
 end

 % marine
 ii=[]; ii=find(ind==-1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y2 y2],'color',clr0(6,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(6,:),'linewidth',2)
 end

 % terrestrial
 ii=[]; ii=find(ind==1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y1 y1],'color',clr0(5,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(5,:),'linewidth',2)
 end
 axis([0 12 -40 0])
 alpha(0.1)

 xlabel('Time (kyr BP)','fontsize',12,'fontweight','normal')
 ylabel('RSL (m)','fontsize',12,'fontweight','normal')
 %title('b. Holocene RSL reconstructions','fontsize',12,'fontweight','normal')

 % do a least squares fit and plot it
 ii=[]; ii=find(ind==0);
 pp=[]; ss=[]; [pp,ss]=polyfit(s(ii),z(ii),1);
 ee=[]; ee=(inv(ss.R)*inv(ss.R)')*((ss.normr)^2)/ss.df;
 % generate solutions
 bb=[];
 for kk=1:1000
  bb(kk,:)=mvnrnd(pp,ee);
 end 
 clear pp ii ee ss
 %X(1,:)=(min(s):10:max(s));
 X(1,:)=0:100:1.2e4;
 X(2,:)=ones(size(X(1,:)));
 c=fill(1e-3*[X(1,:) fliplr(X(1,:))],[prctile(bb*X,97.5) fliplr(prctile(bb*X,2.5))],clr0(8,:)); 
 set(c,'EdgeColor',clr0(8,:),'facealpha',.5,'linewidth',2)
 plot(1e-3*X(1,:),prctile(bb*X,50),'color',clr0(8,:),'linewidth',2)
 legend([{'Index point'};{'Marine limit'};{'Freshwater limit'};{'OLS linear fit'}],'location','south','orientation','vertical'), %legend boxoff

% start making some figures
fig=figure('color','white')
fig.Position(3) = fig.Position(3)*4/3;
fig.Position(4) = fig.Position(4)*2/3;
subplot(1,2,2)

 box on, hold on, grid on
 plot([100 200],[100 200],'color',clr0(7,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(6,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(5,:),'linewidth',2)
 plot([100 200],[100 200],'color',[.5 .5 .5],'linewidth',2)

 ii=[]; ii=find(ind==0);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 c=fill([x1 x2 x2 x1 x1]*1e-3,[y1 y1 y2 y2 y1],'k','linewidth',2);
 set(c,'EdgeColor',clr0(7,:),'facecolor',clr0(7,:))
 end

 % marine
 ii=[]; ii=find(ind==-1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y2 y2],'color',clr0(6,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(6,:),'linewidth',2)
 end

 % terrestrial
 ii=[]; ii=find(ind==1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y1 y1],'color',clr0(5,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(5,:),'linewidth',2)
 end
 axis([0 12 -40 0])
 alpha(0.1)

 xlabel('Time (kyr BP)','fontsize',12,'fontweight','normal')
 ylabel('RSL (m)','fontsize',12,'fontweight','normal')
 %title('b. Holocene RSL reconstructions','fontsize',12,'fontweight','normal')

 [ALPHA1,BETA,GAMMA2]=model1;
 bb=[ALPHA1 BETA];
 X(1,:)=0:100:1.2e4;
 X(2,:)=ones(size(X(1,:)));
 c=fill(1e-3*[X(1,:) fliplr(X(1,:))],[prctile(bb*X,97.5) fliplr(prctile(bb*X,2.5))],clr0(8,:)); 
 set(c,'EdgeColor',clr0(8,:),'facealpha',.5,'linewidth',2)
 plot(1e-3*X(1,:),prctile(bb*X,50),'color',clr0(8,:),'linewidth',2)
 legend([{'Index point'};{'Marine limit'};{'Freshwater limit'};{'Bayes model 1'}],'location','south','orientation','vertical'), %legend boxoff


% start making some figures
fig=figure('color','white')
fig.Position(3) = fig.Position(3)*4/3;
fig.Position(4) = fig.Position(4)*2/3;
subplot(1,2,2)

 box on, hold on, grid on
 plot([100 200],[100 200],'color',clr0(7,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(6,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(5,:),'linewidth',2)
 plot([100 200],[100 200],'color',[.5 .5 .5],'linewidth',2)

 ii=[]; ii=find(ind==0);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 c=fill([x1 x2 x2 x1 x1]*1e-3,[y1 y1 y2 y2 y1],'k','linewidth',2);
 set(c,'EdgeColor',clr0(7,:),'facecolor',clr0(7,:))
 end

 % marine
 ii=[]; ii=find(ind==-1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y2 y2],'color',clr0(6,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(6,:),'linewidth',2)
 end

 % terrestrial
 ii=[]; ii=find(ind==1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y1 y1],'color',clr0(5,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(5,:),'linewidth',2)
 end
 axis([0 12 -40 0])
 alpha(0.1)

 xlabel('Time (kyr BP)','fontsize',12,'fontweight','normal')
 ylabel('RSL (m)','fontsize',12,'fontweight','normal')
 %title('b. Holocene RSL reconstructions','fontsize',12,'fontweight','normal')

 [ALPHA2,BETA,GAMMA2]=model2;
 bb=[ALPHA2 BETA];
 X(1,:)=0:100:1.2e4;
 X(2,:)=ones(size(X(1,:)));
 c=fill(1e-3*[X(1,:) fliplr(X(1,:))],[prctile(bb*X,97.5) fliplr(prctile(bb*X,2.5))],clr0(8,:)); 
 set(c,'EdgeColor',clr0(8,:),'facealpha',.5,'linewidth',2)
 plot(1e-3*X(1,:),prctile(bb*X,50),'color',clr0(8,:),'linewidth',2)
 legend([{'Index point'};{'Marine limit'};{'Freshwater limit'};{'Bayes model 2'}],'location','south','orientation','vertical'), %legend boxoff

% start making some figures
fig=figure('color','white')
fig.Position(3) = fig.Position(3)*4/3;
fig.Position(4) = fig.Position(4)*2/3;
subplot(1,2,2)

 box on, hold on, grid on
 plot([100 200],[100 200],'color',clr0(7,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(6,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(5,:),'linewidth',2)
 plot([100 200],[100 200],'color',[.5 .5 .5],'linewidth',2)

 ii=[]; ii=find(ind==0);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 c=fill([x1 x2 x2 x1 x1]*1e-3,[y1 y1 y2 y2 y1],'k','linewidth',2);
 set(c,'EdgeColor',clr0(7,:),'facecolor',clr0(7,:))
 end

 % marine
 ii=[]; ii=find(ind==-1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y2 y2],'color',clr0(6,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(6,:),'linewidth',2)
 end

 % terrestrial
 ii=[]; ii=find(ind==1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y1 y1],'color',clr0(5,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(5,:),'linewidth',2)
 end
 axis([0 12 -40 0])
 alpha(0.1)

 xlabel('Time (kyr BP)','fontsize',12,'fontweight','normal')
 ylabel('RSL (m)','fontsize',12,'fontweight','normal')
 %title('b. Holocene RSL reconstructions','fontsize',12,'fontweight','normal')

 [Y,TEE,GAMMA,LAMBDA,PHI,YPOST3,YRATE3,tpost]=model3(1);
 c=fill(1e-3*[tpost fliplr(tpost)],[prctile(YPOST3,97.5) fliplr(prctile(YPOST3,2.5))],clr0(8,:)); 
 set(c,'EdgeColor',clr0(8,:),'facealpha',.5,'linewidth',2)
 plot(1e-3*X(1,:),prctile(YPOST3,50),'color',clr0(8,:),'linewidth',2)
 legend([{'Index point'};{'Marine limit'};{'Freshwater limit'};{'Bayes model 3'}],'location','south','orientation','vertical'), %legend boxoff

% start making some figures
fig=figure('color','white')
fig.Position(3) = fig.Position(3)*4/3;
fig.Position(4) = fig.Position(4)*2/3;
subplot(1,2,2)

 box on, hold on, grid on
 plot([100 200],[100 200],'color',clr0(7,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(6,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(5,:),'linewidth',2)
 plot([100 200],[100 200],'color',[.5 .5 .5],'linewidth',2)

 ii=[]; ii=find(ind==0);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 c=fill([x1 x2 x2 x1 x1]*1e-3,[y1 y1 y2 y2 y1],'k','linewidth',2);
 set(c,'EdgeColor',clr0(7,:),'facecolor',clr0(7,:))
 end

 % marine
 ii=[]; ii=find(ind==-1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y2 y2],'color',clr0(6,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(6,:),'linewidth',2)
 end

 % terrestrial
 ii=[]; ii=find(ind==1);
 for kk=1:numel(ii)
 x1=[]; x2=[]; y1=[]; y2=[];
 x1=s(ii(kk))-2*sig(ii(kk));
 x2=s(ii(kk))+2*sig(ii(kk));
 y1=z(ii(kk))-2*zet(ii(kk));
 y2=z(ii(kk))+2*zet(ii(kk));
 plot([x1 x2]*1e-3,[y1 y1],'color',clr0(5,:),'linewidth',2)
 plot(mean([x1 x2])*[1 1]*1e-3,[y1 y2],'color',clr0(5,:),'linewidth',2)
 end
 axis([0 12 -40 0])
 alpha(0.1)

 xlabel('Time (kyr BP)','fontsize',12,'fontweight','normal')
 ylabel('RSL (m)','fontsize',12,'fontweight','normal')
 %title('b. Holocene RSL reconstructions','fontsize',12,'fontweight','normal')

 [Y,TEE,GAMMA,LAMBDA,PHI,YPOST4,YRATE4,tpost]=model4(1);
 c=fill(1e-3*[tpost fliplr(tpost)],[prctile(YPOST4,97.5) fliplr(prctile(YPOST4,2.5))],clr0(8,:)); 
 set(c,'EdgeColor',clr0(8,:),'facealpha',.5,'linewidth',2)
 plot(1e-3*X(1,:),prctile(YPOST4,50),'color',clr0(8,:),'linewidth',2)
 legend([{'Index point'};{'Marine limit'};{'Freshwater limit'};{'Bayes model 4'}],'location','south','orientation','vertical'), %legend boxoff

% start making some figures
fig=figure('color','white')
fig.Position(3) = fig.Position(3)*4/3;
fig.Position(4) = fig.Position(4)*2/3;
subplot(1,2,2)

 box on, hold on, grid on
 plot([100 200],[100 200],'color',clr0(1,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(2,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(3,:),'linewidth',2)
 plot([100 200],[100 200],'color',clr0(4,:),'linewidth',2)


 xlabel('Time (kyr BP)','fontsize',12,'fontweight','normal')
 ylabel('RSL rate (mm/yr)','fontsize',12,'fontweight','normal')
 %title('b. Holocene RSL reconstructions','fontsize',12,'fontweight','normal')

 c=fill(1e-3*[min(tpost) max(tpost) max(tpost) min(tpost) min(tpost)],-1e3*[prctile(ALPHA1,97.5) prctile(ALPHA1,97.5) prctile(ALPHA1,2.5) prctile(ALPHA1,2.5) prctile(ALPHA1,97.5)],clr0(1,:)); 
 set(c,'EdgeColor',clr0(1,:),'facealpha',.8,'linewidth',2)
 plot(1e-3*[min(tpost) max(tpost)],-1e3*median(ALPHA1,1)*[1 1],'color',clr0(1,:),'linewidth',2)

 axis([0 12 -1 5])

 c=fill(1e-3*[min(tpost) max(tpost) max(tpost) min(tpost) min(tpost)],-1e3*[prctile(ALPHA2,97.5) prctile(ALPHA2,97.5) prctile(ALPHA2,2.5) prctile(ALPHA2,2.5) prctile(ALPHA2,97.5)],clr0(2,:)); 
 set(c,'EdgeColor',clr0(2,:),'facealpha',.5,'linewidth',2)
 plot(1e-3*[min(tpost) max(tpost)],-1e3*median(ALPHA2,1)*[1 1],'color',clr0(2,:),'linewidth',2)

 c=fill([tpost fliplr(tpost)]*1e-3,1e3*[prctile(YRATE3,2.5) fliplr(prctile(YRATE3,97.5))],clr0(3,:));
 set(c,'EdgeColor',clr0(3,:),'facealpha',.2,'linewidth',2)
 plot(1e-3*tpost,1e3*median(YRATE3,1),'color',clr0(3,:),'linewidth',2)
 
 c=fill([tpost fliplr(tpost)]*1e-3,1e3*[prctile(YRATE4,2.5) fliplr(prctile(YRATE4,97.5))],clr0(4,:));
 set(c,'EdgeColor',clr0(4,:),'facealpha',.1,'linewidth',2)
 plot(1e-3*tpost,1e3*median(YRATE4,1),'color',clr0(4,:),'linewidth',2)
  legend([{'Model 1'};{'Model 2'};{'Model 3'};{'Model 4'}],'location','northwest','orientation','vertical'), %legend boxoff

