function [ cC, mpsthnom,c_psthnom,all_mrerank,ind_mrerank,crerank,ind_crerank,c_corrm,para,prepostrange,cC_rerank] = cnmfe_yc_mdatalowpassfilter_readmtxt_decimateto1( txtfile,zscoreyn,rawyn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   mrerank: psthnom all videos ranked, index are in ind_mrerank
%   crerank: cells containing matrix, each matrix is the mpsthnom of each
%   video ranked
% c_corrm: cell containing cells, each cell has correlation matrix of
% [wholetrial, beforestim, poststim, delta] for each video


%% initiate txtfile
% txtfile='.txt';
fid=fopen(txtfile);
t=fscanf(fid,'%s',[1]);
expname=t;
t=fscanf(fid,'%s',[1]);

%% initiate variables
mpsthnom=[];
all_mrerank=[];
ind_crerank={};
crerank={};
cC={};
cC_rerank={};
ind=1;
para=struct('ppratio',[],'ppconfi',[]);
c_corrm={};
prepostrange={};
while t
    %% read info
    %     keyboard();
    inmatfile=t
    fprintf(inmatfile)
    t=fscanf(fid,'%s',[1]);
    fr=str2double(t);
    
    t=fscanf(fid,'%s',[1]);
    timestamp=t;
    sind=strfind(timestamp,':');
    minu=str2double(timestamp(1:sind-1));
    seco=str2double(timestamp(sind+1:end));
    timestamp=(minu*60+seco)*fr;
    
    t=fscanf(fid,'%s',[1]);
    pretime=str2double(t)*fr;
    
    t=fscanf(fid,'%s',[1]);
    posttime=str2double(t)*fr;
    
    %% process data
    load(inmatfile,'neuron');
    tpsthnom=[];
%     keyboard();
    if rawyn
        C=cnmfe_lowpassfilter_4hz_20db(neuron.C_raw',fr);
        C=C';
    else
        C=neuron.C;
    end
    temp=C;
    C=[];
    for i=1:size(temp,1)
        C=[C; decimate(temp(i,:),fr)];
    end
    timestamp=timestamp/fr;
    pretime=pretime/fr;
    posttime=posttime/fr;
    fr=1;
    clear temp
    cC{ind}=C;
    if timestamp+posttime<=size(C,2);
        %         keyboard();
        psth=C(:,timestamp-pretime:timestamp+posttime);
%         keyboard();
        for i=1:size(psth,1)
            psth(i,:)=psth(i,:)-mean(psth(i,1:pretime));
            if zscoreyn
%                 keyboard();
                psth(i,:)=psth(i,:)./std(psth(i,1:pretime));
                
            end
        end
%         keyboard();
        
    else % in case the PSTH range is over the size
        psth=C(:,timestamp-pretime:end);
        for i=1:size(psth,1)
            psth(i,:)=psth(i,:)-mean(psth(i,1:pretime));
            if zscoreyn
                psth(i,:)=psth(i,:)./std(psth(i,1:pretime));
            end
        end
        
        psth(:,end:pretime+posttime+1)=NaN;
    end
    
try
    mpsthnom=[mpsthnom;psth]; %matrix with data from all videos
 
catch
    
    keyboard();
end
tpsthnom=[tpsthnom; psth ]; %matrix with data from this video
    
%     keyboard();
    
    [ rerank,ind_rerank] = cnmfe_yc_rankPSTH( tpsthnom,[1:pretime],[pretime:pretime+posttime]);
    %         keyboard();
    crerank{ind}=rerank;
    ind_crerank{ind}=ind_rerank;
    
    C_rerank=C*0;
    for i=1:length(ind_rerank)
        C_rerank(i,:)=C(ind_rerank(i),:);
    end
    cC_rerank{ind}=C_rerank;
    %% get correlation matrix [wholetrial, beforestim, poststim, delta]; consider use GPU
    [allcorr] = analyze_mat2corrmat( rerank );
    [precorr] = analyze_mat2corrmat( rerank(:,1:pretime) );
    [postcorr]=analyze_mat2corrmat( rerank(:,pretime:pretime+posttime) );
    deltacorr=postcorr-precorr;
    temp={allcorr,precorr,postcorr,deltacorr};
    c_corrm{ind}=temp;
    
    
    
    
    
    prepostrange{ind}={[1:pretime],[pretime:pretime+posttime]};
    c_psthnom{ind}=tpsthnom;
    %% re
    
    
    
    ind=ind+1;
    t=fscanf(fid,'%s',[1]);
    
    
end

[ all_mrerank,ind_mrerank ] = cnmfe_yc_rankPSTH( mpsthnom,[1:pretime],[pretime:pretime+posttime]);
    for i2=1:size(all_mrerank,1)
        all_mrerank(i2,:)=all_mrerank(i2,:)-mean(all_mrerank(i2,1:pretime));
    end
prestimF=mpsthnom(:,1:pretime);
poststimF=mpsthnom(:,pretime+1:pretime+posttime);
% keyboard();
para.ppratio=nanmean(poststimF,2)-nanmean(prestimF,2);

para.std=std(mpsthnom')';
para.ppratioz=para.ppratio./para.std;
% keyboard();
end

