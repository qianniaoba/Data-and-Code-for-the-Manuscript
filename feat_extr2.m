function [feat,mpsth]=feat_extr2(spike_data,effcha)
%effcha=[3 8 12 16 17 18 22 24 26];
ncha=size(spike_data,1);
chaid=ismember(1:ncha,effcha);
%cha=cha(chaid);
spikedata=spike_data(chaid,2);
for k=1:size(spikedata,1)
    spiketim=spikedata{k,1};
    t_spk.sk=spiketim;
    feat(k,1)=length(spiketim);
    [spsth(k,:) t]=psth(t_spk,1/20,'n',[0 0.5]);
end
mpsth(:,1)=mean(spsth(k,t>=0 & t<0.5));

