
delete('*PD*')
delete('*.TRB')
delete('*TMP')
delete('*.csv')
delete('*.xlsx')

ballcalib = dir('C:\Users\jcadi\Desktop\CPT\CPT MAT\BALL\*.DAT');
calib = dir('C:\Users\jcadi\Desktop\CPT\CPT MAT\CONE\*.DAT');
numfiles=length(calib);
numfilesball=length(ballcalib);
val = cell(1,numfiles);
for k =1:numfiles
     
    
   
   Aconeproj = .00861113;
Aballproj = .087263022;   diameterball=10.1596238;
   filename = ('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone\calib(k).name');
   r=extractfield(calib,'name');
   z=k;
   v=r(z);
   K=cell2mat(v(:));
   hFig=figure('NumberTitle','off','visible','off','Name',K,'units','normalized','outerposition',[0 0 1 1]);
   h=title(calib(k).name);
   set(h,'Visible','on');
   BASELINES= dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [38 1 38 7]);
 BALLBASELINES= dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\BALL',(ballcalib(k).name)), '%\t', [38 1 38 7]);
 CALFACT0=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [39 2 39 7]);
    CALFACT1=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [40 2 40 7]);
    CALFACT2=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [41 2 41 7]);
    CALFACT3=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [42 2 42 7]);
    CALFACT4=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [43 2 43 7]);
    UNITFACT0=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [44 2 44 7]);
    UNITFACT1=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [45 2 45 7]);
    UNITFACT2=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [46 2 46 7]);
    UNITFACT3=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [47 2 47 7]);
    UNITFACT4=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [48 2 48 7]);
    OFFSET=dlmread(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), '%\t', [49 2 49 7]);
    basesoilmoisture=BASELINES(2);
    baseresistivity=BASELINES(1);
    basetip=BASELINES(4);
    ballbasetip=BALLBASELINES(4);
    basesleeve=BASELINES(5);
    basetemperature=BASELINES(3);
    baseporepressure=BASELINES(6);
    rescalfact0=CALFACT0(1);
    rescalfact1=CALFACT1(1);
    rescalfact2=CALFACT2(1);
    rescalfact3=CALFACT3(1);
    rescalfact4=CALFACT4(1);
    soilmoisturecalfact0=CALFACT0(2);
    soilmoisturecalfact1=CALFACT1(2);
    soilmoisturecalfact2=CALFACT2(2);
    soilmoisturecalfact3=CALFACT3(2);
    soilmoisturecalfact4=CALFACT4(2);
    tempcalfact0=CALFACT0(3);
    tempcalfact1=CALFACT1(3);
    tempcalfact2=CALFACT2(3);
    tempcalfact3=CALFACT3(3);
    tempcalfact4=CALFACT4(3);
    tipcalfact0=CALFACT0(4);
    tipcalfact1=CALFACT1(4);
    tipcalfact2=CALFACT2(4);
    tipcalfact3=CALFACT3(4);
    tipcalfact4=CALFACT4(4);
    sleevecalfact0=CALFACT0(5);
    sleevecalfact1=CALFACT1(5);
    sleevecalfact2=CALFACT2(5);
    sleevecalfact3=CALFACT3(5);
    sleevecalfact4=CALFACT4(5);
    porepressurecalfact0=CALFACT0(6);
    porepressurecalfact1=CALFACT1(6);
    porepressurecalfact2=CALFACT2(6);
    porepressurecalfact3=CALFACT3(6);
    porepressurecalfact4=CALFACT4(6);
    %unit factor assignment
    resunitfact0=UNITFACT0(1);
    resunitfact1=UNITFACT1(1);
    resunitfact2=UNITFACT2(1);
    resunitfact3=UNITFACT3(1);
    resunitfact4=UNITFACT4(1);
    soilmoistunitfact0=UNITFACT0(2);
    soilmoistunitfact1=UNITFACT1(2);
    soilmoistunitfact2=UNITFACT2(2);
    soilmoistunitfact3=UNITFACT3(2);
    soilmoistunitfact4=UNITFACT4(2);
    tempunitfact0=UNITFACT0(3);
    tempunitfact1=UNITFACT1(3);
    tempunitfact2=UNITFACT2(3);
    tempunitfact3=UNITFACT3(3);
    tempunitfact4=UNITFACT4(3);
    tipunitfact0=UNITFACT0(4);
    tipunitfact1=UNITFACT1(4);
    tipunitfact2=UNITFACT2(4);
    tipunitfact3=UNITFACT3(4);
    tipunitfact4=UNITFACT4(4);
    sleeveunitfact0=UNITFACT0(5);
    sleeveunitfact1=UNITFACT1(5);
    sleeveunitfact2=UNITFACT2(5);
    sleeveunitfact3=UNITFACT3(5);
    sleeveunitfact4=UNITFACT4(5);
    porepressureunitfact0=UNITFACT0(6);
    porepressureunitfact1=UNITFACT1(6);
    porepressureunitfact2=UNITFACT2(6);
    porepressureunitfact3=UNITFACT3(6);
    porepressureunitfact4=UNITFACT4(6);
    %channel offset assignment
     resistivityoffset=.275;
    soilmoistureoffset=.275;
    temperatureoffset=OFFSET(3);
    tipoffset=0;
    sleeveoffset=0.1;
    porepressureoffset=0.035;
    
    DATA=readtable(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\Cone',(calib(k).name)), 'HeaderLines',60);
    BALLDATA=readtable(fullfile('C:\Users\jcadi\Desktop\CPT\CPT MAT\BALL',(ballcalib(k).name)), 'HeaderLines',60);
    Depth=[DATA(:,1)];
    BallDepth=[BALLDATA(:,1)];
    
    BallDepth=table2array(BallDepth);
    BallDepth2=BallDepth;
    BallDepthforcorr=BallDepth;
    Depth=table2array(Depth);
    depthend=DATA(end,[1]);
    depthendball=BALLDATA(end,[1]);
    depthendball=table2array(depthendball);
    depthlengthball=length(BallDepth);
    depthlengthball=depthlengthball-1;
    depthend=table2array(depthend);
    depthlength=length(Depth);
    depthlengthcone=depthlength-1;
     TIME=BALLDATA(:,10);
    time=table2array(TIME);
    TIMECONE=DATA(:,10);
    timecone=table2array(TIMECONE);
    speeder = [];
    speedercone = [];
    voltresistivity=DATA(:,2);
    vres=table2array(voltresistivity);
        volttemp=DATA(:,4);
    volttemp=table2array(volttemp);
    volttip=DATA(:,5);
    ballvolttip=BALLDATA(:,5);
    ballvtip=table2array(ballvolttip);
    vtip=table2array(volttip);
     voltsleeve=DATA(:,6);
     vsleeve=table2array(voltsleeve);
    voltpore=DATA(:,7);
     vpore=table2array(voltpore);
     voltmoisture=DATA(:,3);
      vmoist=table2array(voltmoisture);
      moisturecalctemp=(vmoist-soilmoisturecalfact1)/(soilmoisturecalfact2-soilmoisturecalfact1);
     calcrestemp=-2000+(5000./(vres*rescalfact2-rescalfact3));
     res=((calcrestemp.^1.11815E+000)*2.67755E-002);
      temp=(volttemp/(tempcalfact0/1000))-tempcalfact1;
      soilmoist=((moisturecalctemp.^3)*soilmoistunitfact3)+((moisturecalctemp.^2)*soilmoistunitfact2)+(moisturecalctemp*soilmoistunitfact1)+soilmoistunitfact0;
     Tipeng=(vtip-basetip)*(1000/tipcalfact0)*(1/2000)*(1/Aconeproj)*95.76;
    BallTipeng=(ballvtip-ballbasetip)*(1000/0.114064)*(1/2000)*(1/Aconeproj)*(Aconeproj/Aballproj)*95.76;
       % large vane
%     Sleeve=(vsleeve-basesleeve)*(1000/(sleevecalfact0))*((1/0.122718)+(1/(8*0.0278)))*0.04788025889;
       %1.5" small vane
%         Sleeve=(vsleeve-basesleeve)*(1000/(sleevecalfact0))*((1/0.122718)+(1/(8*0.026)))*0.04788025889;

                   Sleeve=(vsleeve-basesleeve)*(1000/(sleevecalfact0))*((1/0.122718)+((0)))*0.04788025889;
% 1 inch sleeve
%  Sleeve=(vsleeve-basesleeve)*(1000/(sleevecalfact0))*((1/0.122718)+(1/(8*0.0208333)))*0.04788025889;
%  2.5 inch
%       Sleeve=(vsleeve-basesleeve)*(1000/(sleevecalfact0))*((1/0.122718)+(1/(8*0.0441)))*0.04788025889;
              
%%0.122718 is area of sleeve---modified to account for sleeve
       %%modifications using fins. one face of fin is 0.0606ft^2, multiply
       %%by 8 to get all faces
     Poreeng=(vpore-baseporepressure)*((1000/porepressurecalfact0))*0.072*95.76;
       DPTHADJRES=Depth-resistivityoffset;
     DPTHADJMOIST=Depth-soilmoistureoffset;
     DPTHADJSLVE=Depth-sleeveoffset;
     DPTHADJPORE=Depth-porepressureoffset;
      for i=1:depthlengthcone
        if i==1
            timeold=0;
            timenew=time(i);
            pushdepthold=0;
            pushdepthnew=Depth(i);
            
        else
             
        
        
        timeold=(time(i))/1000;
        timenew=(time(i+1))/1000;
                pushdepthold=Depth(i);
        pushdepthnew=Depth(i+1);
%         p=i;
%         if timenew == timeold
%             
%             while timenew==timeold && p<depthlengthball
%                 p=p+1
% 
%                 timenew=time(p)/1000;
%                
%             end
%         end
                
        
  
  
        
           


       
        pushspeed=100*((pushdepthnew-pushdepthold)/(timenew-timeold));
        if pushspeed == inf
            pushspeed = nan;
        end
         speeder(i+1,1) = [pushspeed];
  
        
        end
      end
 
%       end
%       for h =1:(length(speeder)-1)
%          if speeder(h+1) == inf
%             speeder(h+1)=speeder(h)
%          end
%       end
%       for h =1:length(speeder)
%           if speeder(h) ==nan
%               BallDepth2(h) = nan
%           end
%       end
%       I = ~isnan(speeder) & ~isnan(BallDepth2);
      
      
%           for i=1:depthlengthcone
%         if i==1
%             timeoldcone=0;
%             timeconenew=timecone(i);
%             pushdeptholdcone=0;
%             pushdepthnew=Depth(i);
%         else
%              
%         
%         
%   
%   
%   
%         timeoldcone=(timecone(i))/1000;
%         timeconenew=(timecone(i+1))/1000;
%         pushdeptholdcone=Depth(i);
%         pushdepthnewcone=Depth(i+1);
%         pushspeedcone=100*((pushdepthnewcone-pushdeptholdcone)/(timeconenew-timeoldcone));
%          speedercone(i+1,1) = [pushspeedcone];
%    
%         end
%     end
%      TOCORR={BallTipeng speeder};
%      BallTipCor=TOCORR{1,1};
%      lenballtimcor=length(BallTipCor);
%      TimeToCor=TOCORR{1,2};
%      lentimetocor=length(TimeToCor);
%      lenspeed=length(speeder)-1;
% 
%      for P=1:lenspeed
%         
%              
%              if speeder(P)<=0 
%                 BallTipCor(P)=nan;
%                 BallDepthforcorr(P)=nan;
%                 BallTipCor(P+1)=nan;
%                 BallDepthforcorr(P+1)=nan;
%                 BallTipCor(P+2)=nan;
%                 BallDepthforcorr(P+2)=nan;
%                
%          
%          
%      
%   
%                 TimeToCor(P)=nan;
%                 TimeToCor(P+1)=nan;
%                 TimeToCor(P+2)=nan;
% 
%                 
%              end
%      end
%          
%      
%              
%      for l=1:lenballtimcor
%          if TimeToCor(l) <=0
% 
%              BallTipCor(l) = nan;
%              BallTipCor(l+1)=nan;
%               BallTipCor(l+2)=nan;
% 
%              TimeToCor(l) = nan;
%              TimeToCor(l+1)=nan;
%              TimeToCor(l+2)=nan;
% 
% %              speeder(l)=nan;
% %              speeder(l+1)=nan;
% %              speeder(l+2)=nan;
% %              speeder(l+3)=nan;
%              BallDepthforcorr(l)=nan;
%              BallDepthforcorr(l+1)=nan;
%              BallDepthforcorr(l+2)=nan;
% 
%                  else
%              
%          end
%      end
%      BallTipCorrected=BallTipCor./((1+(0.11*log((TimeToCor./diameterball)/(2.0/diameterball)))));
% BallTipCor(isnan(BallTipCor))=[];
% BallTipCorrected(isnan(BallTipCorrected))=[];
% BallDepthforcorr(isnan(BallDepthforcorr))=[];
% 
%     
% 
%      PERDIF=-100*(BallTipCor-BallTipCorrected)./BallTipCor;
%    
%      su=BallTipCorrected/10.3;
%    F=fillmissing(su,'movmean',20);
%    D=fillmissing(BallDepthforcorr,'linear');
%    GUNIT=smoothdata(su,'omitnan');
   
%      ax(1)=subplot(1,7,2);                 
%         plot(Poreeng,DPTHADJPORE,'-go',...
%             'Color','black',...
%             'MarkerSize', 3,...
%             'MarkerFaceColor','black');
%        
%         title('Pore Pressure (kPa)')
%         grid on
%         grid minor
%         set(gca,'ydir','reverse')
%         xlim([0,inf]);
%         if depthend > 2.5
%             ylim([0 2])
%             
%         else
%             ylim([0 2])
%             
%         end
%      ax(2) = subplot(1,7,6);
%         plot(BallTipeng,Depth,'-go',...
%             'Color','black',...
%             'MarkerSize', 3,...
%             'MarkerFaceColor','black');
%         title('Tip (kPa)')
%           grid on
%         grid minor
%         set(gca,'ydir','reverse')
%         xlim([0,inf]);
%         if depthend > 2.5
%             ylim([0 2])
%             
%         else
%             ylim([0 2])
%             
%         end
%        
% ax(3)= subplot(1,7,3);
%         plot(Sleeve,DPTHADJSLVE,'-go',...
%             'Color','black',...
%             'MarkerSize', 3,...
%             'MarkerFaceColor','black');
%         title ('fs (kPa)');
%           grid on
%         grid minor
%         set(gca,'ydir','reverse')
%         xlim([0,inf]);
%         if depthend > 2.5
%             ylim([0 2])
%             
%         else
%             ylim([0 2])
%             
%         end
%  ax4= subplot(1,7,4);
%         
%         plot(soilmoist, DPTHADJMOIST, '-go',...
%             'Color','black',...
%             'MarkerSize', 3,...
%             'MarkerFaceColor','black');
%            grid on
%         grid minor
%         set(gca,'ydir','reverse')
% %         xlim([inf,100]);
%         title('Soil Moisture');
%         if depthend > 2.5
%             ylim([0 2]);
%             
%         else
%             ylim([0 2]);
%             
%         end
%   ax5=subplot(1,7,5);
%         
%         plot(res, DPTHADJRES, '-go',...
%             'Color','black',...
%             'MarkerSize', 3,...
%             'MarkerFaceColor','black');
%              grid on
%         grid minor
%         set(gca,'ydir','reverse')
%         xlim([0,100]);
%         title('Resistivity');
%         xlim([0 5])
%         if depthend > 2.5
%             ylim([0 2]);
%             
%         else
%             ylim([0 2]);
%             
%         end
%          ax5=subplot(1,7,1);
%         
%         plot(speeder,Depth,  '-go',...
%             'Color','black',...
%             'MarkerSize', 3,...
%             'MarkerFaceColor','black');
%              grid on
%         grid minor
%         set(gca,'ydir','reverse')
%         xlim([-inf,inf]);
%         title('Push Speed (cm/s)');
%         if depthend > 2.5
%             ylim([0 2])
%             
%         else
%             ylim([0 2])
%             
%         end
% %         ax5=subplot(1,7,7);
% %         
% %         plot(speedercone, Depth,  '-go',...
% %             'Color','black',...
% %             'MarkerSize', 3,...
% %             'MarkerFaceColor','black');
% %              grid on
% %         grid minor
% %         set(gca,'ydir','reverse')
% %         xlim([0,10]);
% %         title('Push Speed Ball(cm/s)');
% %         if depthend > 2.5
% %             ylim([0 4.5])
% %             
% %         else
% %             ylim([0 2.5])
% %             
% %         end
%         suptitle(K);
        
        
clear NEW;
clear vq1;
clear x_idx;
clear idx;
clear xx;
clear depth;
clear AA;
clear BB;
clear CC;
clear DD;
clear EE;
clear BF;
clear bq1;
clear pq1;
clear mq1;
clear rq1;
AA=[Depth BallTipeng];
BB=[Depth Tipeng];
CC=[Depth Poreeng];
DD=[Depth res];
EE=[Depth soilmoist];
FF=[Depth Sleeve];
BF=[Depth temp];
A=rmmissing(AA);
depth=AA(:,1);
% depth=table2array(depth);
ballre=AA(:,2);
porere=CC(:,2);
resre=DD(:,2);
moisre=EE(:,2);
sleevere=FF(:,2);
tipre=BB(:,2);
tempre=BF(:,2);
% sleeve=table2array(sleeve);
[Z,ia,idx] = unique(depth,'last');
val = accumarray(idx,ballre,[],@max);
val2 = accumarray(idx,porere,[],@max);
val3 = accumarray(idx,resre,[],@max);
val4 = accumarray(idx,moisre,[],@max);
val5 = accumarray(idx,sleevere,[],@max);
val6 = accumarray(idx,tipre,[],@max);
val7 = accumarray(idx,tempre,[],@max);
% for f=2:length(depth);
%     if depth(f)==depth(f-1)
%         depth(f)=depth(f)+0.000001
%     else
%     end
% end
% AA=[depth Sleeve];
% [~,idx]=unique(AA(:,1));
% 
% out=AA(idx,:)
xx=[0:0.005:depth(end)];
xx=xx';
% depth=out(:,1);
% depth=table2array(depth);
% sleeve=out(:,2);
% sleeve=table2array(sleeve);
%%linearly re-interpolate at 0.5cm intervals

vq1=interp1(Z,val,xx,'linear');
vq2=interp1(Z,val2,xx,'linear');
vq3=interp1(Z,val3,xx,'linear');
vq4=interp1(Z,val4,xx,'linear');
vq5=interp1(Z,val5,xx,'linear');
vq6=interp1(Z,val6,xx,'linear');
vq7=interp1(Z,val7,xx,'linear');
NEW=[xx,vq1];
NEW2=[xx,vq2];
depthinterpolated=xx;
depthinterpolatedmoist=depthinterpolated-0.27;
ballinterpolated=vq1;
poreinterpolated=vq2;
resisinterpolated=vq3;
moistureinterpolated=vq4;
sleeveinterpolated=vq5;
tipinterpolated=vq6;
tempinterpolated=vq7;
ballsmoothed=smoothdata(ballinterpolated,'sgolay');
poresmoothed=smoothdata(poreinterpolated,'sgolay');
ressmoothed=smoothdata(resisinterpolated,'sgolay');
moistsmoothed=smoothdata(moistureinterpolated,'sgolay');
sleevesmoothed=smoothdata(sleeveinterpolated,'sgolay');
tipsmoothed=smoothdata(tipinterpolated,'sgolay');
temperaturesmoothed=smoothdata(tempinterpolated,'sgolay');
      header={'Depth(m)','Push Speed (cm/s)','Fins-Sleeve(kPa)','Pore Pressure(kPa)','Soil Moisture(%)','Resistivity(Ohm-meters)','Temperature (C)','Ball Tip(UNC)(kPa)','Interpolated Depth(m)','Interpolated Ball Tip(kPa)','Smoothed Ball','Pore Pressure Smoothed','Depth Interp. For SMR','Resistivity Smoothed','Soil Moisture Smoothed','Sleeve Smoothed','Conical Tip Smoothed','temperature smoothed','conicalunprocessed'};
      A=Depth;
      B=speeder;
      C=Sleeve;
      D=Poreeng;
      E=soilmoist;
      F=res;
      GA=temp;
      G=BallTipeng;
      H=depthinterpolated;
      I=ballinterpolated;
      J=ballsmoothed;
      J2=poresmoothed;
      J3=depthinterpolatedmoist;
      J4=ressmoothed;
      J5=moistsmoothed;
      J6=sleevesmoothed;
      J7=tipsmoothed;
      J8=temperaturesmoothed;
      JH=Tipeng;
%       I=speeder;
      M=padcat(A,B,C,D,E,F,GA,G,H,I,J,J2,J3,J4,J5,J6,J7,J8,JH);
      xlswrite('brian.xlsx',header,calib(k).name);
      xlswrite('brian.xlsx',M,calib(k).name,'A2');
%       xlswritefig(hFig,'MAR 2019 CPT.xlsx',calib(k).name,'T1')

end
