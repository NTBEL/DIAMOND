clear
clc

addpath( genpath( 'DIAMOND DIR' ) );

% Code run time for N= 3000 is about 56.255795 seconds.

N=3000; % # of nanobubble signals for analysis

% initialization
signalresultamp=zeros(N,1);
signalresultauc=zeros(N,1);
signalresultlf=zeros(N,1);
signalresultprob=zeros(1,1);
signalresultSNR_ratio=zeros(N,1);

bubble_lf=zeros(1,N);
bubble_amp=zeros(1,N);
bubble_auc=zeros(1,N);
bubble_prob=zeros(1,N);
bubble_lf_down=zeros(1,N);
bubble_lf_up=zeros(1,N);
SNR_ratio=zeros(1,N);
Taillevel=zeros(1,N);

% read the data files
startN =0;
for n=(startN+1:startN+N) 
     if n-1<=9
         filename=['C2Trace0000' num2str(n-1) '.txt'];
     elseif n-1<=99
         filename=['C2Trace000' num2str(n-1) '.txt']; 
     elseif n-1<=999
         filename=['C2Trace00' num2str(n-1) '.txt'];
     else 
         filename=['C2Trace0' num2str(n-1) '.txt'];
     end

format long;
    [fileID,mm]=fopen(filename,'r');
    mydata=textscan(fileID,'%f,%f', 'HeaderLines', 5);
    rawdata= [mydata{1,1},mydata{1,2}];
    fclose(fileID);
format long;
deltatime= rawdata(2,1)-rawdata(1,1);
rawdata(:,2)=rawdata(:,2)-mean(rawdata(50:500,2));
time=rawdata(:,1);

% signal filering
windowSize=20;
b=(1/windowSize)*ones(1,windowSize);
a=1e-4;
Y=filter(b,a,rawdata(:,2));

% calculate the signal to noise ratio
SNR_ratio(1,n)=abs(min(Y(500:5000)))/abs(3*std(Y(4800:5000)));
Taillevel(1,n)=mean(Y(4900:5002));

% signal to noise ratio-based signal classification
if SNR_ratio(1,n)<15
    if SNR_ratio(1,n)<3.7
        bubble_lf(1,n)=0;
        bubble_amp(1,n)=0;
        bubble_auc(1,n)=0;
        bubble_prob(1,n)=0;
        minvalue=0;
        positiondown=1;
        positionup=1; 
        time(positiondown)=0;
        time(positionup)=0;
        Y(positiondown:positionup)=0;
    else
    Y=detrend(Y);
    Y=Y-mean(Y(50:400));
    bubble_prob(1,n)=1;
    [minvalue,position]=min(Y);
    positiondown=position-1;
    positionup=position+1;    
    while Y(positiondown-1)<=mean(Y(50:500))&& positiondown-1>500
        positiondown=positiondown-1;
    end    
    while Y(positionup+1)<=mean(Y(4000:5000))&& positionup+1<5000
        positionup=positionup+1;
    end
    end
else
    [minvalue,position]=min(Y);
    bubble_prob(1,n)=1;
    positiondown=500;% nanobubble signal starting point
    positionup=5002;% index should be greater than end point
    while Y(positiondown+1)>(mean(Y(50:500))-3*abs(std(Y(50:500))))
        positiondown=positiondown+1;
    end
    if Taillevel(1,n)<(min(Y)-mean(Y(50:500)))*0.5
       positionup=length(Y);
    else
        k=(mean(Y(4000:5000))-min(Y(4000:5000)))/std(Y(4000:5000));
        k=ceil(abs(k));
        while Y(positionup-1)>(mean(Y(4900:5000))-2.5*k*abs(std(Y(4500:5000))))
            positionup=positionup-1;
        end
    end
end
bubble_lf_down(1,n)=time(positiondown)*10^9;
bubble_lf_up(1,n)=time(positionup)*10^9;
bubble_lf(1,n)=(time(positionup)-time(positiondown))*10^9;
bubble_amp(1,n)=abs(minvalue);
bubble_auc(1,n)=sum(abs(Y(positiondown:positionup)))*deltatime*10^9;
end

signalresultlf_down=(bubble_lf_down)';
signalresultlf_up=(bubble_lf_up)';
signalresultlf(:,1)=(bubble_lf)';
signalresultamp(:,1)=(bubble_amp)';
signalresultauc(:,1)=(bubble_auc)';
signalresultSNR_ratio(:,1)=(SNR_ratio)';
signalresultprob(1,1)=sum(bubble_prob)/N*100;

csvwrite('2. Target_amplitude.csv',signalresultamp);
csvwrite('2. Target_area under the curve.csv',signalresultauc);

f=figure;
f.Position = [500 300 350 350];
scatter(signalresultamp,signalresultauc,'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
              'MarkerFaceColor',[0.9290 0.6940 0.1250],...
              'LineWidth',1.5)
xlim([(min(signalresultamp)-50) (max(signalresultamp)+50)])
xlabel('Amplitude (mV)')
ylabel('AUC (mV.ns)')
legend('Target sample')
set(gca,'FontSize',12)