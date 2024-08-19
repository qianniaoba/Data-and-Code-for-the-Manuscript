clc
clear
close all



load('IPPG_signal_sample.mat')
load('Contact_based_device_results.mat')

HR_range1 = 3*60;
HR_range2 = 15*60;


s=IPPG_signal_sample;
Res_C=Contact_based_device_results;
  
m=size(s,2);
B=fft(s);
C=abs(B);

Freq = 1:m;
Freq = (Freq-1)/m*frA;  
mask =( Freq >(HR_range1/60) & Freq <(HR_range2/60));     

B(~mask)=0;
Q=real(ifft(B));

s=Q;

D=abs(fft(s));
% fs = 30;
fs = frA;
t = (1/fs):(1/fs):(length(s)/fs);

%% CMOR Wavelet Transform

wavename='cmor80-6';

totalscal=256*32;     

Fc=centfrq(wavename); 


c=2*Fc*totalscal;
scals=c./(1:totalscal);
f=scal2frq(scals,wavename,1/fs); 
coefs=cwt(s,scals,wavename); 
W=abs(coefs);


W(1:(floor(HR_range1 / (frA / 2 / totalscal * 60))),:)=0;
W((ceil(HR_range2 / (frA / 2 / totalscal * 60))+1):end,:)=0;
Copy_W=W;
i=0;


h=figure('InvertHardcopy','off','Color',[1 1 1]);
currentPosition = get(h, 'Position');
width = 800; 
height = 300; 
set(h, 'Position', [currentPosition(1), currentPosition(2), width, height]);
axes1 = axes;
hold(axes1,'on');

% imagesc(t,f*60,abs(coefs));
imagesc(t,f*60,W);
% axis([min(t) max(t) min(f) 4*60]);

axis([min(t) max(t) min(f*60) max(f*60)]);

% axis([min(t) max(t) min(f*60) 240]);

colormap jet;
set(gca,'YDir','normal')

xlim(axes1,[0 size(IPPG_signal_sample,2)/frA]);
% xticks(0:600:size(IPPG_signal_sample,2)/frA);
set(gca, 'XTick',[0:600:size(IPPG_signal_sample,2)/frA]);

ylim(axes1,[HR_range1 HR_range2]);


xlabel('Time (s)','FontSize',12);

ylabel('Heart rate (bpm)','FontSize',12);

box(axes1,'on');
set(axes1,'FontSize',12);

%% Heart Rate Extraction

while i~=m
    
    [row,line]=find(Copy_W==max(max(Copy_W)));

    range = ceil(10 / (frA / 2 / totalscal * 60));

    temp1=row;
    temp2=row;
    RR=zeros(1,m);
    PKS=zeros(1,m);
    RR(line) = f(temp1);
    PKS(line)=max(max(Copy_W));
    Indexmax=max(max(Copy_W));
    for i=line-1:-1:1
       [pks,locs]= findpeaks(W(temp1-range:temp1+range,i),'MinPeakDistance',2*range-1);
       temp1 = temp1 - range + locs -1;
       RR(i) = f(temp1);
       PKS(i)=pks;
%        if (pks/Indexmax)<0.4
%            Copy_W(:,i:2*line-i)=0;
%           
%            break
%        end
    end
    if i==1
        for i=line+1:m
           [pks,locs]= findpeaks(W(temp2-range:temp2+range,i),'MinPeakDistance',2*range-1);
           temp2 = temp2 - range + locs -1;
           RR(i) = f(temp2);
           PKS(i)=pks;
        end
    end
end

% figure
% % imagesc(t,f*60,abs(coefs));
% imagesc(t,f*60,W);
% axis([min(t) max(t) min(f) 4*60]);
% colormap jet;
% set(gca,'YDir','normal')
% 
% colorbar;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Energy spectrum plot');


% ft=fittype(['fourier','3']);
% opts=fitoptions(ft);
% [fitres,gof]=fit(t',HR',ft,opts);
% yf=fitres(t);
% figure
% % plot(t,yf,'r','LineWidth',3)
% % hold on;

% RR_HR=RR*60;
% plot(t,RR_HR,'r','LineWidth',2)
% % ylim(axes1,[200 500]);
% axis([-inf,inf,200,500])
% hold on;

% figure
% 
% RR_HR=RR*60;
% plot(t,RR_HR,'r','LineWidth',2)
% 
% axis([-inf,inf,200,500])
% hold on;
% plot(t,Res_C,'b','LineWidth',2)
% xlabel('Time (s)');
% ylabel('Heart rate (bpm)');

h=figure('InvertHardcopy','off','Color',[1 1 1]);
currentPosition = get(h, 'Position');
width = 800; 
height = 300; 
set(h, 'Position', [currentPosition(1), currentPosition(2), width, height]);
axes1 = axes;
hold(axes1,'on');

HR_RR=RR*60;
plot(t,HR_RR,'r','LineWidth',2)
hold on;
plot(t,Res_C,'b','LineWidth',2)
%  ylim(axes1,[200 500]);
% axis([-inf,inf,200,500])
hold on;
xlim(axes1,[0 size(IPPG_signal_sample,2)/frA]);
% xticks(0:600:size(IPPG_signal_sample,2)/frA);
set(gca, 'XTick',[0:600:size(IPPG_signal_sample,2)/frA]);


ylim(axes1,[200 500]);


xlabel('Time (s)','FontSize',12);

ylabel('Heart rate (bpm)','FontSize',12);

box(axes1,'on');
set(axes1,'FontSize',12);
%%
% X=HR_RR(0.5*frA/2:0.5*frA:end);
% Y=Res_C(0.5*frA/2:0.5*frA:end);
% 
% ME=mean(abs(X-Y))
% 
% SD=std((X-Y))
% 
% RMSE=rmse(X,Y)
% 
% r=corrcoef(X,Y)






