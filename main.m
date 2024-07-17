clc
clear
close all

load('video_image.mat')

k = 0 ;

yuzhi = 0.3;
frA = 60; % frames per second (fps)

HR_range1 = 3*60;
HR_range2 = 15*60;
    
R_G=cell(1,size(video_image,2));

for i=1:size(video_image,2)
    k = k+1
    out=video_image{1, i};
    figure(1);imshow(out);
    % figure(1);imshow(out(:,:,3)<yuzhi);

    rows=size(out,1);
    columns=size(out,2); 
   
    rank=[];
    n=0;
     for i=1:rows
        for j=1:columns
            
            n = n +1 ;
            rank(1,n)=out(i,j,3);           

            rank(2,n)=i;
            rank(3,n)=j;
            rank(4,n)=out(i,j,2);
                  
        end
     end
     rank=rank';
     rank=sortrows(rank,1);
     rank=rank';
   

     count=sum(rank(1,:)<yuzhi);
     R_G(1,k)={sort(rank(4,1:count))};

     % temp_rank=rank(:,1:count);

     % temp_rank=temp_rank';
     % temp_rank=sortrows(temp_rank,4);
     % temp_rank=temp_rank';

     % if sum(abs(temp_rank(4,:)-sort(rank(4,1:count))))~=0
     %     fprintf('---')
     % end

end

%% Set Sequence Interval 
temp_P=[];
for j=1:size(R_G,2)
    temp_P=[temp_P size(R_G{1,j},2)];
end

RR_Y=NaN(size(R_G,2),max(temp_P) );

for j=1:size(R_G,2)
    RR_Y(j,1:size(R_G{1,j},2))=R_G{1,j};
end

C_C=nanmean(RR_Y);

[xData, yData] = prepareCurveData( [], C_C );


ft = fittype( 'fourier5' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';


[fitresult, gof] = fit( xData, yData, ft, opts );

% figure;plot(diff(fitresult(xData),1));hold on;plot([1 size(xData,1)],[min(diff(fitresult(xData),1)) min(diff(fitresult(xData),1))]*3)
NiHe_diff_1=diff(fitresult(xData),1);

% WZ=find(NiHe_diff_1==min(NiHe_diff_1)); 

WZ=find(NiHe_diff_1==min(NiHe_diff_1)); 

for i=WZ:-1:1
    if NiHe_diff_1(i)>(min(NiHe_diff_1)*3)
        range_1=i+1;
        break

    end
end


for i=WZ:1:size(NiHe_diff_1,1)
    if NiHe_diff_1(i)>(min(NiHe_diff_1)*3)
        range_2=i-1;
        break

    end
end

% figure;plot(fitresult(xData));hold on;plot(range_1,fitresult(range_1),'.',MarkerSize=15,Color=[1 0 0]);hold on;plot(range_2,fitresult(range_2),'.',MarkerSize=15,Color=[1 0 0])


h=figure('InvertHardcopy','off','Color',[1 1 1]);
currentPosition = get(h, 'Position');
width = 550; 
height = 400; 
set(h, 'Position', [currentPosition(1), currentPosition(2), width, height]);
axes1 = axes;
hold(axes1,'on');

for i=1:size(fitresult(xData),1)
    hold on
    
    plot(i,fitresult(i),'.','Color',[0 fitresult(i) 0]);
end
hold on;plot([range_1:range_2],fitresult([range_1:range_2]),'.','Color',[248,149,136]/255);
hold on;plot(range_1,fitresult(range_1),'.','MarkerSize',15,'Color',[1 0 0]);
hold on;plot(range_2,fitresult(range_2),'.','MarkerSize',15,'Color',[1 0 0]);
ylim(axes1,[0 0.4]);
xlabel('Cumulative pixel number','FontSize',12);

ylabel('Pixel mean value','FontSize',12);

box(axes1,'on');
set(axes1,'FontSize',12);

h=figure('InvertHardcopy','off','Color',[1 1 1]);
currentPosition = get(h, 'Position');
width = 800; 
height = 300; 
set(h, 'Position', [currentPosition(1), currentPosition(2), width, height]);
axes1 = axes;
hold(axes1,'on');

    temp_C=unique(RR_Y(:));
    temp_R=temp_C(temp_C>=0);
    for i=1:size(temp_R,1)
        [a b]=find(RR_Y==temp_R(i,1));
        % for j=1:size(a,1)
            hold on
            plot(a,ones(1,size(a,1))*temp_R(i,1),'.','Color',[0 temp_R(i,1) 0])

        % end

    end

for i=1:1:size(RR_Y,1)
    % i
    hold on

    hold on
    plot(ones(1,range_2-range_1+1)*i,RR_Y(i,range_1:range_2),'.','Color',[248,149,136]/255)

    hold on
    plot(i,RR_Y(i,range_1),'.','Color',[1 0 0])
    
    hold on
    plot(i,RR_Y(i,range_2),'.','Color',[1 0 0])
   

end

xlim(axes1,[0 size(video_image,2)]);

xlabel('Video frame','FontSize',12);
ylabel('Pixel value','FontSize',12);

box(axes1,'on');
set(axes1,'FontSize',12);

%% Raw IPPG Signal
h=figure('InvertHardcopy','off','Color',[1 1 1]);
currentPosition = get(h, 'Position');
width = 800; 
height = 300; 
set(h, 'Position', [currentPosition(1), currentPosition(2), width, height]);
axes1 = axes;
hold(axes1,'on');

plot(mean(RR_Y(:,range_1:range_2)'),'LineWidth',1,'Color',[1 0 0])

xlim(axes1,[0 size(video_image,2)]);

ylim(axes1,[0.130 0.175]);
xlabel('Video frame','FontSize',12);

ylabel('Pixel mean value','FontSize',12);

box(axes1,'on');
set(axes1,'FontSize',12);

%%

% s=mean(RR_Y(:,range_1:range_2)');
% 
% 
% 
% 
% %%    
% 
% m=size(s,2);
% B=fft(s);
% C=abs(B);
% 
% 
% Freq = 1:m;
% 
% Freq = (Freq-1)/m*frA;  
% 
% mask =( Freq >(180/60) & Freq <(900/60));   
% %
% % mask =( Freq >(20/60) & Freq <(900/60));   
% 
% B(~mask)=0;
% Q=real(ifft(B));
% 
% s=Q;
% 
% D=abs(fft(s));
% figure;plot([1:size(s,2)]/60,s)



% fs = 30;
%% CMOR Wavelet Transform
fs = frA;

s=filt(mean(RR_Y(:,range_1:range_2)'),60,3,15); 
t = (1/fs):(1/fs):(length(s)/fs);
m=size(s,2);

wavename='cmor5-6';

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

xlim(axes1,[0 size(video_image,2)/frA]);

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
% plot(t,yf,'r','LineWidth',3)
% hold on;

h=figure('InvertHardcopy','off','Color',[1 1 1]);
currentPosition = get(h, 'Position');
width = 800; 
height = 300; 
set(h, 'Position', [currentPosition(1), currentPosition(2), width, height]);
axes1 = axes;
hold(axes1,'on');

HR_RR=RR*60;
plot(t,HR_RR,'r','LineWidth',3)
%  ylim(axes1,[200 500]);
% axis([-inf,inf,200,500])
hold on;
xlim(axes1,[0 size(video_image,2)/frA]);

ylim(axes1,[200 500]);


xlabel('Time (s)','FontSize',12);

ylabel('Heart rate (bpm)','FontSize',12);

box(axes1,'on');
set(axes1,'FontSize',12);




