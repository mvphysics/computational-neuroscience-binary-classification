

% Calling eeglab and closing its window
eeglab; close;
% Setting parameters
clear all
zalpha=1.96;

% Loading data
[b1,b1s]=sload('B0701T.gdf');
[b2,b2s]=sload('B0702T.gdf');
% Resampling from 250Hz to 100Hz
b1=resample(b1,100,250);
b2=resample(b2,100,250);

% Creating event windows
b1_type=b1s.EVENT.TYP;
b2_type=b2s.EVENT.TYP;
b1_pos=b1s.EVENT.POS;
b2_pos=b2s.EVENT.POS;
check=0;
for i=1:length(b1_type)
    if b1_type(i)==768
        if b1_type(i+1)==769 || b1_type(i+1)==770
            new_pos=floor(b1_pos(i+1)/2.5);
            new_row=[b1_type(i+1),new_pos,new_pos+399];
            if check==1
                b1_data=[b1_data;new_row];
            else
                b1_data=new_row;
                check=1;
            end
        end
    end
end
clear b1_type b1_pos
check=0;
for i=1:length(b2_type)
    if b2_type(i)==768
        if b2_type(i+1)==769 || b2_type(i+1)==770
            new_pos=floor(b2_pos(i+1)/2.5);
            new_row=[b2_type(i+1),new_pos,new_pos+399];
            if check==1
                b2_data=[b2_data;new_row];
            else
                b2_data=new_row;
                check=1;
            end
        end
    end
end
clear b2_type b2_pos new_row new_pos

% Plotting sample signals
l_idx=find(b1_data(1:end,1)==769,1);
r_idx=find(b1_data(1:end,1)==770,1);
color=[0.8,0.4,0.85];
fig_1=figure(1);
t = tiledlayout(3,2,'TileSpacing','Compact');
subplot(3,2,1);
plot(b1(b1_data(l_idx,2):b1_data(l_idx,3),1),'Color',color);
title('Cue onset LEFT | C3');
xlabel('100Hz*sec');
xlim([0 400]);
subplot(3,2,3);
plot(b1(b1_data(l_idx,2):b1_data(l_idx,3),2),'Color',color);
title('Cue onset LEFT | Cz');
xlabel('100Hz*sec');
xlim([0 400]);
subplot(3,2,5);
plot(b1(b1_data(l_idx,2):b1_data(l_idx,3),3),'Color',color);
title('Cue onset LEFT | C4')
xlim([0 400]);
xlabel('100Hz*sec');
subplot(3,2,2);
plot(b1(b1_data(r_idx,2):b1_data(r_idx,3),1),'Color',color);
title('Cue onset RIGHT | C3');
xlabel('100Hz*sec');
xlim([0 400]);
subplot(3,2,4);
plot(b1(b1_data(r_idx,2):b1_data(r_idx,3),2),'Color',color);
title('Cue onset RIGHT | Cz');
xlabel('100Hz*sec');
xlim([0 400]);
subplot(3,2,6);
plot(b1(b1_data(r_idx,2):b1_data(r_idx,3),3),'Color',color);
title('Cue onset RIGHT | C4');
xlabel('100Hz*sec');
xlim([0 400]);
sgtitle('First Session');

l_idx=find(b2_data(1:end,1)==769,1);
r_idx=find(b2_data(1:end,1)==770,1);
color=[0.8,0.4,0.85];
fig_2=figure(2);
t = tiledlayout(3,2,'TileSpacing','Compact');
subplot(3,2,1);
plot(b2(b2_data(l_idx,2):b2_data(l_idx,3),1),'Color',color);
title('Cue onset LEFT | C3');
xlabel('100Hz*sec');
xlim([0 400]);
subplot(3,2,3);
plot(b2(b2_data(l_idx,2):b2_data(l_idx,3),2),'Color',color);
title('Cue onset LEFT | Cz');
xlabel('100Hz*sec');
xlim([0 400]);
subplot(3,2,5);
plot(b2(b2_data(l_idx,2):b2_data(l_idx,3),3),'Color',color);
title('Cue onset LEFT | C4')
xlim([0 400]);
xlabel('100Hz*sec');
subplot(3,2,2);
plot(b2(b2_data(r_idx,2):b2_data(r_idx,3),1),'Color',color);
title('Cue onset RIGHT | C3');
xlabel('100Hz*sec');
xlim([0 400]);
subplot(3,2,4);
plot(b2(b2_data(r_idx,2):b2_data(r_idx,3),2),'Color',color);
title('Cue onset RIGHT | Cz');
xlabel('100Hz*sec');
xlim([0 400]);
subplot(3,2,6);
plot(b2(b2_data(r_idx,2):b2_data(r_idx,3),3),'Color',color);
title('Cue onset RIGHT | C4');
xlabel('100Hz*sec');
xlim([0 400]);
sgtitle('Second Session');


% Calculating metrics for single channel processing
% First for the left data
n_lag=50;
l_idx=find(b1_data(1:end,1)==769);
C3_check=0;
Cz_check=0;
C4_check=0;
for i=1:length(l_idx)
    [acf,lag]=autocorr(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),1),NumLags=n_lag);
    mean_C3=mean(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),1));
    std_C3=std(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),1));
    sk_C3=skewness(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),1));
    kur_C3=kurtosis(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),1));
    iqr_C3=iqr(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),1));
    [act_C3,mob_C3,compl_C3]=hjorth(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),1));
    if C3_check==1
        l_autocorr_C3=[l_autocorr_C3,acf];
    else
        l_autocorr_C3=[lag,acf];
        C3_check=1;
    end
    [acf,lag]=autocorr(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),2),NumLags=n_lag);
    mean_Cz=mean(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),2));
    std_Cz=std(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),2));
    sk_Cz=skewness(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),2));
    kur_Cz=kurtosis(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),2));
    iqr_Cz=iqr(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),2));
    [act_Cz,mob_Cz,compl_Cz]=hjorth(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),2));
    if Cz_check==1
        l_autocorr_Cz=[l_autocorr_Cz,acf];
    else
        l_autocorr_Cz=[lag,acf];
        Cz_check=1;
    end
    [acf,lag]=autocorr(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),3),NumLags=n_lag);
    mean_C4=mean(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),3));
    std_C4=std(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),3));
    sk_C4=skewness(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),3));
    kur_C4=kurtosis(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),3));
    iqr_C4=iqr(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),3));
    [act_C4,mob_C4,compl_C4]=hjorth(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),3));
    if C4_check==1
        l_autocorr_C4=[l_autocorr_C4,acf];
        new_row=[mean_C3,std_C3,sk_C3,kur_C3,iqr_C3,act_C3,mob_C3,compl_C3,mean_Cz,std_Cz,sk_Cz,kur_Cz,iqr_Cz,act_Cz,mob_Cz,compl_Cz,mean_C4,std_C4,sk_C4,kur_C4,iqr_C4,act_C4,mob_C4,compl_C4];
        left_data=[left_data;new_row];
    else
        l_autocorr_C4=[lag,acf];
        left_data=[mean_C3,std_C3,sk_C3,kur_C3,iqr_C3,act_C3,mob_C3,compl_C3,mean_Cz,std_Cz,sk_Cz,kur_Cz,iqr_Cz,act_Cz,mob_Cz,compl_Cz,mean_C4,std_C4,sk_C4,kur_C4,iqr_C4,act_C4,mob_C4,compl_C4];
        C4_check=1;
    end
end
l_idx=find(b2_data(1:end,1)==769);
for i=1:length(l_idx)
    [acf,lag]=autocorr(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),1),NumLags=n_lag);
    mean_C3=mean(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),1));
    std_C3=std(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),1));
    sk_C3=skewness(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),1));
    kur_C3=kurtosis(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),1));
    iqr_C3=iqr(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),1));
    [act_C3,mob_C3,compl_C3]=hjorth(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),1));
    l_autocorr_C3=[l_autocorr_C3,acf];
    [acf,lag]=autocorr(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),2),NumLags=n_lag);
    mean_Cz=mean(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),2));
    std_Cz=std(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),2));
    sk_Cz=skewness(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),2));
    kur_Cz=kurtosis(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),2));
    iqr_Cz=iqr(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),2));
    [act_Cz,mob_Cz,compl_Cz]=hjorth(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),2));
    l_autocorr_Cz=[l_autocorr_Cz,acf];
    [acf,lag]=autocorr(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),3),NumLags=n_lag);
    mean_C4=mean(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),3));
    std_C4=std(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),3));
    sk_C4=skewness(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),3));
    kur_C4=kurtosis(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),3));
    iqr_C4=iqr(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),3));
    [act_C4,mob_C4,compl_C4]=hjorth(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),3));
    l_autocorr_C4=[l_autocorr_C4,acf];
    new_row=[mean_C3,std_C3,sk_C3,kur_C3,iqr_C3,act_C3,mob_C3,compl_C3,mean_Cz,std_Cz,sk_Cz,kur_Cz,iqr_Cz,act_Cz,mob_Cz,compl_Cz,mean_C4,std_C4,sk_C4,kur_C4,iqr_C4,act_C4,mob_C4,compl_C4];
    left_data=[left_data;new_row];
end

% Then for the right data
r_idx=find(b1_data(1:end,1)==770);
C3_check=0;
Cz_check=0;
C4_check=0;
for i=1:length(r_idx)
    [acf,lag]=autocorr(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),1),NumLags=n_lag);
    mean_C3=mean(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),1));
    std_C3=std(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),1));
    sk_C3=skewness(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),1));
    kur_C3=kurtosis(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),1));
    iqr_C3=iqr(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),1));
    [act_C3,mob_C3,compl_C3]=hjorth(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),1));
    if C3_check==1
        r_autocorr_C3=[r_autocorr_C3,acf];
    else
        r_autocorr_C3=[lag,acf];
        C3_check=1;
    end
    [acf,lag]=autocorr(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),2),NumLags=n_lag);
    mean_Cz=mean(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),2));
    std_Cz=std(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),2));
    sk_Cz=skewness(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),2));
    kur_Cz=kurtosis(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),2));
    iqr_Cz=iqr(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),2));
    [act_Cz,mob_Cz,compl_Cz]=hjorth(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),2));
    if Cz_check==1
        r_autocorr_Cz=[r_autocorr_Cz,acf];
    else
        r_autocorr_Cz=[lag,acf];
        Cz_check=1;
    end
    [acf,lag]=autocorr(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),3),NumLags=n_lag);
    mean_C4=mean(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),3));
    std_C4=std(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),3));
    sk_C4=skewness(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),3));
    kur_C4=kurtosis(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),3));
    iqr_C4=iqr(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),3));
    [act_C4,mob_C4,compl_C4]=hjorth(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),3));
    if C4_check==1
        r_autocorr_C4=[r_autocorr_C4,acf];
        new_row=[mean_C3,std_C3,sk_C3,kur_C3,iqr_C3,act_C3,mob_C3,compl_C3,mean_Cz,std_Cz,sk_Cz,kur_Cz,iqr_Cz,act_Cz,mob_Cz,compl_Cz,mean_C4,std_C4,sk_C4,kur_C4,iqr_C4,act_C4,mob_C4,compl_C4];
        right_data=[right_data;new_row];
    else
        r_autocorr_C4=[lag,acf];
        right_data=[mean_C3,std_C3,sk_C3,kur_C3,iqr_C3,act_C3,mob_C3,compl_C3,mean_Cz,std_Cz,sk_Cz,kur_Cz,iqr_Cz,act_Cz,mob_Cz,compl_Cz,mean_C4,std_C4,sk_C4,kur_C4,iqr_C4,act_C4,mob_C4,compl_C4];
        C4_check=1;
    end
end
r_idx=find(b2_data(1:end,1)==770);
for i=1:length(r_idx)
    [acf,lag]=autocorr(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),1),NumLags=n_lag);
    mean_C3=mean(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),1));
    std_C3=std(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),1));
    sk_C3=skewness(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),1));
    kur_C3=kurtosis(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),1));
    iqr_C3=iqr(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),1));
    [act_C3,mob_C3,compl_C3]=hjorth(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),1));
    r_autocorr_C3=[r_autocorr_C3,acf];
    [acf,lag]=autocorr(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),2),NumLags=n_lag);
    mean_Cz=mean(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),2));
    std_Cz=std(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),2));
    sk_Cz=skewness(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),2));
    kur_Cz=kurtosis(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),2));
    iqr_Cz=iqr(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),2));
    [act_Cz,mob_Cz,compl_Cz]=hjorth(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),2));
    r_autocorr_Cz=[r_autocorr_Cz,acf];
    [acf,lag]=autocorr(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),3),NumLags=n_lag);
    mean_C4=mean(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),3));
    std_C4=std(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),3));
    sk_C4=skewness(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),3));
    kur_C4=kurtosis(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),3));
    iqr_C4=iqr(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),3));
    [act_C4,mob_C4,compl_C4]=hjorth(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),3));
    r_autocorr_C4=[r_autocorr_C4,acf];
    new_row=[mean_C3,std_C3,sk_C3,kur_C3,iqr_C3,act_C3,mob_C3,compl_C3,mean_Cz,std_Cz,sk_Cz,kur_Cz,iqr_Cz,act_Cz,mob_Cz,compl_Cz,mean_C4,std_C4,sk_C4,kur_C4,iqr_C4,act_C4,mob_C4,compl_C4];
    right_data=[right_data;new_row];
end

clear l_idx r_idx acf mean_C3 std_C3 sk_C3 kur_C3 iqr_C3 act_C3 mob_C3 compl_C3 mean_Cz std_Cz sk_Cz kur_Cz iqr_Cz act_Cz mob_Cz compl_Cz mean_C4 std_C4 sk_C4 kur_C4 iqr_C4 act_C4 mob_C4 compl_C4

% Converting arrays to tables for easier access
left_data=array2table(left_data,'VariableNames',{'C3_MEAN','C3_STD','C3_SKEWNESS','C3_KURTOSIS','C3_INTERQUARTILE','C3_HJORTH_ACT','C3_HJORTH_MOB','C3_HJORTH_COMPLEX', ...
    'Cz_MEAN','Cz_STD','Cz_SKEWNESS','Cz_KURTOSIS','Cz_INTERQUARTILE','Cz_HJORTH_ACT','Cz_HJORTH_MOB','Cz_HJORTH_COMPLEX', ...
    'C4_MEAN','C4_STD','C4_SKEWNESS','C4_KURTOSIS','C4_INTERQUARTILE','C4_HJORTH_ACT','C4_HJORTH_MOB','C4_HJORTH_COMPLEX',});

right_data=array2table(right_data,'VariableNames',{'C3_MEAN','C3_STD','C3_SKEWNESS','C3_KURTOSIS','C3_INTERQUARTILE','C3_HJORTH_ACT','C3_HJORTH_MOB','C3_HJORTH_COMPLEX', ...
    'Cz_MEAN','Cz_STD','Cz_SKEWNESS','Cz_KURTOSIS','Cz_INTERQUARTILE','Cz_HJORTH_ACT','Cz_HJORTH_MOB','Cz_HJORTH_COMPLEX', ...
    'C4_MEAN','C4_STD','C4_SKEWNESS','C4_KURTOSIS','C4_INTERQUARTILE','C4_HJORTH_ACT','C4_HJORTH_MOB','C4_HJORTH_COMPLEX',});


% Plotting data
fig_3=figure(3);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data),left_data.C3_MEAN,'Color','magenta');
hold on
plot(1:height(left_data),left_data.C3_STD,'Color','blue');
plot(1:height(left_data),left_data.C3_SKEWNESS,'Color','green');
plot(1:height(left_data),left_data.C3_KURTOSIS,'Color','yellow');
plot(1:height(left_data),left_data.C3_INTERQUARTILE,'Color','red');
plot(1:height(left_data),(mean(left_data.C3_MEAN)*ones(height(left_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C3_STD)*ones(height(left_data),1)),'Color','blue','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C3_SKEWNESS)*ones(height(left_data),1)),'Color','green','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C3_KURTOSIS)*ones(height(left_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C3_INTERQUARTILE)*ones(height(left_data),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data)])
ylim([-2 12])
title('Cue onset LEFT')
legend(left_data.Properties.VariableNames{1:5},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data),right_data.C3_MEAN,'Color','magenta');
hold on
plot(1:height(right_data),right_data.C3_STD,'Color','blue');
plot(1:height(right_data),right_data.C3_SKEWNESS,'Color','green');
plot(1:height(right_data),right_data.C3_KURTOSIS,'Color','yellow');
plot(1:height(right_data),right_data.C3_INTERQUARTILE,'Color','red');
plot(1:height(right_data),(mean(right_data.C3_MEAN)*ones(height(right_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C3_STD)*ones(height(right_data),1)),'Color','blue','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C3_SKEWNESS)*ones(height(right_data),1)),'Color','green','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C3_KURTOSIS)*ones(height(right_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C3_INTERQUARTILE)*ones(height(right_data),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(right_data)])
ylim([-2 12])
title('Cue onset RIGHT')
legend(right_data.Properties.VariableNames{1:5},'Interpreter', 'none')
hold off

fig_4=figure(4);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data),left_data.C3_HJORTH_ACT,'Color','magenta');
hold on
plot(1:height(left_data),left_data.C3_HJORTH_MOB,'Color','cyan');
plot(1:height(left_data),left_data.C3_HJORTH_COMPLEX,'Color','blue');
plot(1:height(left_data),(mean(left_data.C3_HJORTH_ACT)*ones(height(left_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C3_HJORTH_MOB)*ones(height(left_data),1)),'Color','cyan','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C3_HJORTH_COMPLEX)*ones(height(left_data),1)),'Color','blue','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data)])
ylim([0 40])
title('Cue onset LEFT')
legend(left_data.Properties.VariableNames{6:8},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data),right_data.C3_HJORTH_ACT,'Color','magenta');
hold on
plot(1:height(right_data),right_data.C3_HJORTH_MOB,'Color','cyan');
plot(1:height(right_data),right_data.C3_HJORTH_COMPLEX,'Color','blue');
plot(1:height(right_data),(mean(right_data.C3_HJORTH_ACT)*ones(height(right_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C3_HJORTH_MOB)*ones(height(right_data),1)),'Color','cyan','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C3_HJORTH_COMPLEX)*ones(height(right_data),1)),'Color','blue','LineWidth',2);
xlabel('Sample index')
xlim([1 height(right_data)])
ylim([0 40])
title('Cue onset RIGHT')
legend(right_data.Properties.VariableNames{6:8},'Interpreter', 'none')
hold off

fig_5=figure(5);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data),left_data.Cz_MEAN,'Color','magenta');
hold on
plot(1:height(left_data),left_data.Cz_STD,'Color','blue');
plot(1:height(left_data),left_data.Cz_SKEWNESS,'Color','green');
plot(1:height(left_data),left_data.Cz_KURTOSIS,'Color','yellow');
plot(1:height(left_data),left_data.Cz_INTERQUARTILE,'Color','red');
plot(1:height(left_data),(mean(left_data.Cz_MEAN)*ones(height(left_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data),(mean(left_data.Cz_STD)*ones(height(left_data),1)),'Color','blue','LineWidth',2);
plot(1:height(left_data),(mean(left_data.Cz_SKEWNESS)*ones(height(left_data),1)),'Color','green','LineWidth',2);
plot(1:height(left_data),(mean(left_data.Cz_KURTOSIS)*ones(height(left_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(left_data),(mean(left_data.Cz_INTERQUARTILE)*ones(height(left_data),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data)])
ylim([-2 12])
title('Cue onset LEFT')
legend(left_data.Properties.VariableNames{9:13},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data),right_data.Cz_MEAN,'Color','magenta');
hold on
plot(1:height(right_data),right_data.Cz_STD,'Color','blue');
plot(1:height(right_data),right_data.Cz_SKEWNESS,'Color','green');
plot(1:height(right_data),right_data.Cz_KURTOSIS,'Color','yellow');
plot(1:height(right_data),right_data.Cz_INTERQUARTILE,'Color','red');
plot(1:height(right_data),(mean(right_data.Cz_MEAN)*ones(height(right_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data),(mean(right_data.Cz_STD)*ones(height(right_data),1)),'Color','blue','LineWidth',2);
plot(1:height(right_data),(mean(right_data.Cz_SKEWNESS)*ones(height(right_data),1)),'Color','green','LineWidth',2);
plot(1:height(right_data),(mean(right_data.Cz_KURTOSIS)*ones(height(right_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(right_data),(mean(right_data.Cz_INTERQUARTILE)*ones(height(right_data),1)),'Color','red','LineWidth',2);

xlabel('Sample index')
xlim([1 height(right_data)])
ylim([-2 12])
title('Cue onset RIGHT')
legend(right_data.Properties.VariableNames{9:13},'Interpreter', 'none')
hold off

fig_6=figure(6);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data),left_data.Cz_HJORTH_ACT,'Color','magenta');
hold on
plot(1:height(left_data),left_data.Cz_HJORTH_MOB,'Color','cyan');
plot(1:height(left_data),left_data.Cz_HJORTH_COMPLEX,'Color','blue');
plot(1:height(left_data),(mean(left_data.Cz_HJORTH_ACT)*ones(height(left_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data),(mean(left_data.Cz_HJORTH_MOB)*ones(height(left_data),1)),'Color','cyan','LineWidth',2);
plot(1:height(left_data),(mean(left_data.Cz_HJORTH_COMPLEX)*ones(height(left_data),1)),'Color','blue','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data)])
ylim([0 40])
title('Cue onset LEFT')
legend(left_data.Properties.VariableNames{14:16},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data),right_data.Cz_HJORTH_ACT,'Color','magenta');
hold on
plot(1:height(right_data),right_data.Cz_HJORTH_MOB,'Color','cyan');
plot(1:height(right_data),right_data.Cz_HJORTH_COMPLEX,'Color','blue');
plot(1:height(right_data),(mean(right_data.Cz_HJORTH_ACT)*ones(height(right_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data),(mean(right_data.Cz_HJORTH_MOB)*ones(height(right_data),1)),'Color','cyan','LineWidth',2);
plot(1:height(right_data),(mean(right_data.Cz_HJORTH_COMPLEX)*ones(height(right_data),1)),'Color','blue','LineWidth',2);
xlabel('Sample index')
xlim([1 height(right_data)])
ylim([0 40])
title('Cue onset RIGHT')
legend(right_data.Properties.VariableNames{14:16},'Interpreter', 'none')
hold off

fig_7=figure(7);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data),left_data.C4_MEAN,'Color','magenta');
hold on
plot(1:height(left_data),left_data.C4_STD,'Color','blue');
plot(1:height(left_data),left_data.C4_SKEWNESS,'Color','green');
plot(1:height(left_data),left_data.C4_KURTOSIS,'Color','yellow');
plot(1:height(left_data),left_data.C4_INTERQUARTILE,'Color','red');
plot(1:height(left_data),(mean(left_data.C4_MEAN)*ones(height(left_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4_STD)*ones(height(left_data),1)),'Color','blue','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4_SKEWNESS)*ones(height(left_data),1)),'Color','green','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4_KURTOSIS)*ones(height(left_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4_INTERQUARTILE)*ones(height(left_data),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data)])
ylim([-2 12])
title('Cue onset LEFT')
legend(left_data.Properties.VariableNames{17:22},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data),right_data.C4_MEAN,'Color','magenta');
hold on
plot(1:height(right_data),right_data.C4_STD,'Color','blue');
plot(1:height(right_data),right_data.C4_SKEWNESS,'Color','green');
plot(1:height(right_data),right_data.C4_KURTOSIS,'Color','yellow');
plot(1:height(right_data),right_data.C4_INTERQUARTILE,'Color','red');
plot(1:height(right_data),(mean(right_data.C4_MEAN)*ones(height(right_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4_STD)*ones(height(right_data),1)),'Color','blue','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4_SKEWNESS)*ones(height(right_data),1)),'Color','green','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4_KURTOSIS)*ones(height(right_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4_INTERQUARTILE)*ones(height(right_data),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(right_data)])
ylim([-2 12])
title('Cue onset RIGHT')
legend(right_data.Properties.VariableNames{17:22},'Interpreter', 'none')
hold off

fig_8=figure(8);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data),left_data.C4_HJORTH_ACT,'Color','magenta');
hold on
plot(1:height(left_data),left_data.C4_HJORTH_MOB,'Color','cyan');
plot(1:height(left_data),left_data.C4_HJORTH_COMPLEX,'Color','blue');
plot(1:height(left_data),(mean(left_data.C4_HJORTH_ACT)*ones(height(left_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4_HJORTH_MOB)*ones(height(left_data),1)),'Color','cyan','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4_HJORTH_COMPLEX)*ones(height(left_data),1)),'Color','blue','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data)])
ylim([0 40])
title('Cue onset LEFT')
legend(left_data.Properties.VariableNames{22:24},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data),right_data.C4_HJORTH_ACT,'Color','magenta');
hold on
plot(1:height(right_data),right_data.C4_HJORTH_MOB,'Color','cyan');
plot(1:height(right_data),right_data.C4_HJORTH_COMPLEX,'Color','blue');
plot(1:height(right_data),(mean(right_data.C4_HJORTH_ACT)*ones(height(right_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4_HJORTH_MOB)*ones(height(right_data),1)),'Color','cyan','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4_HJORTH_COMPLEX)*ones(height(right_data),1)),'Color','blue','LineWidth',2);
xlabel('Sample index')
xlim([1 height(right_data)])
ylim([0 40])
title('Cue onset RIGHT')
legend(right_data.Properties.VariableNames{22:24},'Interpreter', 'none')
hold off

% Calculating mean and standard deviation for lags
C3_check=0;
Cz_check=0;
C4_check=0;
for i=2:n_lag+1
    m=mean(l_autocorr_C3(i,2:end));
    s=std(l_autocorr_C3(i,2:end));
    new_row=[m,s];
    if C3_check==1
        l_m_std_C3=[l_m_std_C3;new_row];
    else
        l_m_std_C3=new_row;
        C3_check=1;
    end
    m=mean(l_autocorr_Cz(i,2:end));
    s=std(l_autocorr_Cz(i,2:end));
    new_row=[m,s];
    if Cz_check==1
        l_m_std_Cz=[l_m_std_Cz;new_row];
    else
        l_m_std_Cz=new_row;
        Cz_check=1;
    end
    m=mean(l_autocorr_C3(i,2:end));
    s=std(l_autocorr_C4(i,2:end));
    new_row=[m,s];
    if C4_check==1
        l_m_std_C4=[l_m_std_C4;new_row];
    else
        l_m_std_C4=new_row;
        C4_check=1;
    end
end
C3_check=0;
Cz_check=0;
C4_check=0;
for i=2:n_lag+1
    m=mean(r_autocorr_C3(i,2:end));
    s=std(r_autocorr_C3(i,2:end));
    new_row=[m,s];
    if C3_check==1
        r_m_std_C3=[r_m_std_C3;new_row];
    else
        r_m_std_C3=new_row;
        C3_check=1;
    end
    m=mean(r_autocorr_Cz(i,2:end));
    s=std(r_autocorr_Cz(i,2:end));
    new_row=[m,s];
    if Cz_check==1
        r_m_std_Cz=[r_m_std_Cz;new_row];
    else
        r_m_std_Cz=new_row;
        Cz_check=1;
    end
    m=mean(r_autocorr_C3(i,2:end));
    s=std(r_autocorr_C4(i,2:end));
    new_row=[m,s];
    if C4_check==1
        r_m_std_C4=[r_m_std_C4;new_row];
    else
        r_m_std_C4=new_row;
        C4_check=1;
    end
end

clear C_3check C4_check Cz_check  m s new_row

dif_C3=l_m_std_C3(1:end,1)-r_m_std_C3(1:end,1);
dif_Cz=l_m_std_Cz(1:end,1)-r_m_std_Cz(1:end,1);
dif_C4=l_m_std_C4(1:end,1)-r_m_std_C4(1:end,1);
fig_9=figure(9);
t = tiledlayout(3,2,'TileSpacing','Compact');
subplot(3,2,1);
plot(lag(2:end),l_m_std_C3(1:end,1),'Color','red')
hold on
plot(lag(2:end),r_m_std_C3(1:end,1),'Color','green')
xlabel('Lag')
title('C3 channel | MEAN')
legend('Cue onset LEFT','Cue onset RIGTH')
xlim([1 length(lag)-1])
hold off
subplot(3,2,2);
plot(lag(2:end),dif_C3,'Color','blue')
hold on
plot(lag(2:end),r_m_std_C3(1:end,2),'Color','green','LineStyle','--')
plot(lag(2:end),l_m_std_C3(1:end,2),'Color','red','LineStyle','--')
xlabel('Lag')
title('C3 channel | STD & DIFFERENCE')
legend('Mean difference','Cue onset LEFT','Cue onset RIGTH')
xlim([1 length(lag)-1])
hold off
subplot(3,2,3);
plot(lag(2:end),l_m_std_Cz(1:end,1),'Color','red')
hold on
plot(lag(2:end),r_m_std_Cz(1:end,1),'Color','green')
xlabel('Lag')
title('Cz channel | MEAN')
legend('Cue onset LEFT','Cue onset RIGTH')
xlim([1 length(lag)-1])
hold off
subplot(3,2,4);
plot(lag(2:end),dif_Cz,'Color','blue')
hold on
plot(lag(2:end),l_m_std_Cz(1:end,2),'Color','red','LineStyle','--')
plot(lag(2:end),r_m_std_Cz(1:end,2),'Color','green','LineStyle','--')
xlabel('Lag')
title('Cz channel | STD & DIFFERENCE')
legend('Mean difference','Cue onset LEFT','Cue onset RIGTH')
xlim([1 length(lag)-1])
hold off
subplot(3,2,5);
plot(lag(2:end),l_m_std_C4(1:end,1),'Color','red')
hold on
plot(lag(2:end),r_m_std_C4(1:end,1),'Color','green')
xlabel('Lag')
title('C4 channel | MEAN')
legend('Cue onset LEFT','Cue onset RIGTH')
xlim([1 length(lag)-1])
hold off
subplot(3,2,6);
plot(lag(2:end),dif_C4,'Color','blue')
hold on
plot(lag(2:end),l_m_std_C4(1:end,2),'Color','red','LineStyle','--')
plot(lag(2:end),r_m_std_C4(1:end,2),'Color','green','LineStyle','--')
xlabel('Lag')
title('C4 channel | STD & DIFFERENCE')
legend('Mean difference','Cue onset LEFT','Cue onset RIGTH')
xlim([1 length(lag)-1])
hold off

% Multichannel analysis
% White noise process and calculation of cross correlation Granger causality and conditional Granger causality
p=3;
check=0;
gci_lag=15;
l_idx=find(b1_data(1:end,1)==769);
for i=1:length(l_idx)
    C3=fitAR(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),1),p);
    Cz=fitAR(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),2),p);
    C4=fitAR(b1(b1_data(l_idx(i),2):b1_data(l_idx(i),3),3),p);
    if i==1
        fig_10=figure(10);
        t = tiledlayout(3,2,'TileSpacing','Compact');
        subplot(3,2,1);
        plot(l_autocorr_C3(1:end,1),l_autocorr_C3(1:end,2),'Color',[0.5,0.1,0.5])
        hold on
        plot(l_autocorr_C3(1:end,1),((zalpha/sqrt(400))*ones(length(l_autocorr_C3(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        plot(l_autocorr_C3(1:end,1),(-(zalpha/sqrt(400))*ones(length(l_autocorr_C3(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        title('Autocorrelation of raw data | C3')
        xlabel('lag')
        hold off
        [acf,lag]=autocorr(C3,NumLags=n_lag);
        subplot(3,2,2);
        plot(lag,acf,'Color',[0.5,0.1,0.5])
        hold on
        plot(lag,((zalpha/sqrt(400))*ones(length(l_autocorr_C3(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        plot(lag,(-(zalpha/sqrt(400))*ones(length(l_autocorr_C3(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        title('Autocorrelation after white noise process (p=3) | C3')
        xlabel('lag')
        subplot(3,2,3);
        plot(l_autocorr_Cz(1:end,1),l_autocorr_Cz(1:end,2),'Color',[0.5,0.1,0.5])
        hold on
        plot(l_autocorr_Cz(1:end,1),((zalpha/sqrt(400))*ones(length(l_autocorr_Cz(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        plot(l_autocorr_Cz(1:end,1),(-(zalpha/sqrt(400))*ones(length(l_autocorr_Cz(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        title('Autocorrelation of raw data | Cz')
        xlabel('lag')
        hold off
        [acf,lag]=autocorr(Cz,NumLags=n_lag);
        subplot(3,2,4);
        plot(lag,acf,'Color',[0.5,0.1,0.5])
        hold on
        plot(lag,((zalpha/sqrt(400))*ones(length(l_autocorr_Cz(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        plot(lag,(-(zalpha/sqrt(400))*ones(length(l_autocorr_Cz(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        title('Autocorrelation after white noise process (p=3) | Cz')
        xlabel('lag')
        hold off
        subplot(3,2,5);
        [acf,lag]=autocorr(C4,NumLags=n_lag);
        plot(l_autocorr_C4(1:end,1),l_autocorr_C4(1:end,2),'Color',[0.5,0.1,0.5])
        hold on
        plot(l_autocorr_C4(1:end,1),((zalpha/sqrt(400))*ones(length(l_autocorr_C4(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        plot(l_autocorr_C4(1:end,1),(-(zalpha/sqrt(400))*ones(length(l_autocorr_C4(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        title('Autocorrelation of raw data | C4')
        xlabel('lag')
        hold off
        [acf,lag]=autocorr(C4,NumLags=n_lag);
        subplot(3,2,6);
        plot(lag,acf,'Color',[0.5,0.1,0.5])
        hold on
        plot(lag,((zalpha/sqrt(400))*ones(length(l_autocorr_C4(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        plot(lag,(-(zalpha/sqrt(400))*ones(length(l_autocorr_C4(1:end,1)),1)),'c--','Color',[0.9,0.5,0.9])
        title('Autocorrelation after white noise process (p=3) | C4')
        xlabel('lag')
        hold off
    end
    cor=corrcoef([C3,Cz,C4]);
    if check==1
        new_row=[cor(1,2),cor(2,3),cor(1,3)];
        cr_cor=[cr_cor;new_row];
    else
        cr_cor=[cor(1,2),cor(2,3),cor(1,3)];
    end
    cor=GCI([C3,Cz,C4],gci_lag);
    if check==1
        new_row=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
        gci_cor=[gci_cor;new_row];
    else
        gci_cor=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
    end
    cor=CGCI([C3,Cz,C4],gci_lag);
    if check==1
        new_row=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
        cgci_cor=[cgci_cor;new_row];
    else
        cgci_cor=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
        check=1;
    end
end

l_idx=find(b2_data(1:end,1)==769);
for i=1:length(l_idx)
    C3=fitAR(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),1),p);
    Cz=fitAR(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),2),p);
    C4=fitAR(b2(b2_data(l_idx(i),2):b2_data(l_idx(i),3),3),p);
    cor=corrcoef([C3,Cz,C4]);
    new_row=[cor(1,2),cor(2,3),cor(1,3)];
    cr_cor=[cr_cor;new_row];
    cor=GCI([C3,Cz,C4],gci_lag);
    new_row=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
    gci_cor=[gci_cor;new_row];
    cor=CGCI([C3,Cz,C4],gci_lag);
    new_row=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
    cgci_cor=[cgci_cor;new_row];
end

cr_cor=array2table(cr_cor,'VariableNames',{'C3Cz_Cross_Corr','CzC4_Cross_Corr','C3C4_Cross_Corr'});
gci_cor=array2table(gci_cor,'VariableNames',{'C3Cz_GCI','C3C4_GCI','CzC3_GCI','CzC4_GCI','C4C3_GCI','C4Cz_GCI'});
cgci_cor=array2table(cgci_cor,'VariableNames',{'C3Cz_CGCI','C3C4_CGCI','CzC3_CGCI','CzC4_CGCI','C4C3_CGCI','C4Cz_CGCI'});
left_data=[left_data,cr_cor,gci_cor,cgci_cor];
clear cr_cor cgi_cor cgci_cor

% Adding autocorrelations to the data
C=l_autocorr_C3.';
C=C(2:end,2:11);
C=array2table(C,'VariableNames',{'C3_Lag1','C3_Lag2','C3_Lag3','C3_Lag4','C3_Lag5','C3_Lag6','C3_Lag7','C3_Lag8','C3_Lag9','C3_Lag10'});
left_data=[left_data,C];
C=l_autocorr_Cz.';
C=C(2:end,2:11);
C=array2table(C,'VariableNames',{'Cz_Lag1','Cz_Lag2','Cz_Lag3','Cz_Lag4','Cz_Lag5','Cz_Lag6','Cz_Lag7','Cz_Lag8','Cz_Lag9','Cz_Lag10'});
left_data=[left_data,C];
C=l_autocorr_C4.';
C=C(2:end,2:11);
C=array2table(C,'VariableNames',{'C4_Lag1','C4_Lag2','C4_Lag3','C4_Lag4','C4_Lag5','C4_Lag6','C4_Lag7','C4_Lag8','C4_Lag9','C4_Lag10'});
left_data=[left_data,C];


check=0;
r_idx=find(b1_data(1:end,1)==770);
for i=1:length(r_idx)
    C3=fitAR(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),1),p);
    Cz=fitAR(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),2),p);
    C4=fitAR(b1(b1_data(r_idx(i),2):b1_data(r_idx(i),3),3),p);
    cor=corrcoef([C3,Cz,C4]);
    if check==1
        new_row=[cor(1,2),cor(2,3),cor(1,3)];
        cr_cor=[cr_cor;new_row];
    else
        cr_cor=[cor(1,2),cor(2,3),cor(1,3)];
    end
    cor=GCI([C3,Cz,C4],gci_lag);
    if check==1
        new_row=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
        gci_cor=[gci_cor;new_row];
    else
        gci_cor=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
    end
    cor=CGCI([C3,Cz,C4],gci_lag);
    if check==1
        new_row=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
        cgci_cor=[cgci_cor;new_row];
    else
        cgci_cor=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
        check=1;
    end
end

r_idx=find(b2_data(1:end,1)==770);
for i=1:length(r_idx)
    C3=fitAR(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),1),p);
    Cz=fitAR(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),2),p);
    C4=fitAR(b2(b2_data(r_idx(i),2):b2_data(r_idx(i),3),3),p);
    cor=corrcoef([C3,Cz,C4]);
    new_row=[cor(1,2),cor(2,3),cor(1,3)];
    cr_cor=[cr_cor;new_row];
    cor=GCI([C3,Cz,C4],gci_lag);
    new_row=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
    gci_cor=[gci_cor;new_row];
    cor=CGCI([C3,Cz,C4],gci_lag);
    new_row=[cor(1,2),cor(1,3),cor(2,1),cor(2,3),cor(3,1),cor(3,2)];
    cgci_cor=[cgci_cor;new_row];
end

cr_cor=array2table(cr_cor,'VariableNames',{'C3Cz_Cross_Corr','CzC4_Cross_Corr','C3C4_Cross_Corr'});
gci_cor=array2table(gci_cor,'VariableNames',{'C3Cz_GCI','C3C4_GCI','CzC3_GCI','CzC4_GCI','C4C3_GCI','C4Cz_GCI'});
cgci_cor=array2table(cgci_cor,'VariableNames',{'C3Cz_CGCI','C3C4_CGCI','CzC3_CGCI','CzC4_CGCI','C4C3_CGCI','C4Cz_CGCI'});
right_data=[right_data,cr_cor,gci_cor,cgci_cor];
clear cr_cor cgi_cor cgci_cor

% Adding autocorrelations to the data
C=r_autocorr_C3.';
C=C(2:end,2:11);
C=array2table(C,'VariableNames',{'C3_Lag1','C3_Lag2','C3_Lag3','C3_Lag4','C3_Lag5','C3_Lag6','C3_Lag7','C3_Lag8','C3_Lag9','C3_Lag10'});
right_data=[right_data,C];
C=r_autocorr_Cz.';
C=C(2:end,2:11);
C=array2table(C,'VariableNames',{'Cz_Lag1','Cz_Lag2','Cz_Lag3','Cz_Lag4','Cz_Lag5','Cz_Lag6','Cz_Lag7','Cz_Lag8','Cz_Lag9','Cz_Lag10'});
right_data=[right_data,C];
C=r_autocorr_C4.';
C=C(2:end,2:11);
C=array2table(C,'VariableNames',{'C4_Lag1','C4_Lag2','C4_Lag3','C4_Lag4','C4_Lag5','C4_Lag6','C4_Lag7','C4_Lag8','C4_Lag9','C4_Lag10'});
right_data=[right_data,C];

% Plotting multichannel features
orange=[1,0.7,0];
fig_11=figure(11);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data),left_data.C3Cz_Cross_Corr,'Color','magenta');
hold on
plot(1:height(left_data),left_data.CzC4_Cross_Corr,'Color','cyan');
plot(1:height(left_data),left_data.C3C4_Cross_Corr,'Color','blue');
plot(1:height(left_data),(mean(left_data.C3Cz_Cross_Corr)*ones(height(left_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data),(mean(left_data.CzC4_Cross_Corr)*ones(height(left_data),1)),'Color','cyan','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C3C4_Cross_Corr)*ones(height(left_data),1)),'Color','blue','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data)])
ylim([0 0.6])
title('Cue onset LEFT')
legend(left_data.Properties.VariableNames{25:27},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data),right_data.C3Cz_Cross_Corr,'Color','magenta');
hold on
plot(1:height(right_data),right_data.CzC4_Cross_Corr,'Color','cyan');
plot(1:height(right_data),right_data.C3C4_Cross_Corr,'Color','blue');
plot(1:height(right_data),(mean(right_data.C3Cz_Cross_Corr)*ones(height(right_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data),(mean(right_data.CzC4_Cross_Corr)*ones(height(right_data),1)),'Color','cyan','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C3C4_Cross_Corr)*ones(height(right_data),1)),'Color','blue','LineWidth',2);
xlabel('Sample index')
xlim([1 height(right_data)])
ylim([0 0.6])
title('Cue onset RIGHT')
legend(right_data.Properties.VariableNames{25:27},'Interpreter', 'none')
sgtitle('Cross Correlation')
hold off

fig_12=figure(12);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data),left_data.C3Cz_GCI,'Color','magenta');
hold on
plot(1:height(left_data),left_data.C3C4_GCI,'Color','blue');
plot(1:height(left_data),left_data.CzC3_GCI,'Color','green');
plot(1:height(left_data),left_data.CzC4_GCI,'Color','yellow');
plot(1:height(left_data),left_data.C4C3_GCI,'Color',orange);
plot(1:height(left_data),left_data.C4Cz_GCI,'Color','red');
plot(1:height(left_data),(mean(left_data.C3Cz_GCI)*ones(height(left_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C3C4_GCI)*ones(height(left_data),1)),'Color','blue','LineWidth',2);
plot(1:height(left_data),(mean(left_data.CzC3_GCI)*ones(height(left_data),1)),'Color','green','LineWidth',2);
plot(1:height(left_data),(mean(left_data.CzC4_GCI)*ones(height(left_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4C3_GCI)*ones(height(left_data),1)),'Color',orange,'LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4Cz_GCI)*ones(height(left_data),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data)])
ylim([0 0.2])
title('Cue onset LEFT')
legend(left_data.Properties.VariableNames{28:33},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data),right_data.C3Cz_GCI,'Color','magenta');
hold on
plot(1:height(right_data),right_data.C3C4_GCI,'Color','blue');
plot(1:height(right_data),right_data.CzC3_GCI,'Color','green');
plot(1:height(right_data),right_data.CzC4_GCI,'Color','yellow');
plot(1:height(right_data),right_data.C4C3_GCI,'Color',orange);
plot(1:height(right_data),right_data.C4Cz_GCI,'Color','red');
plot(1:height(right_data),(mean(right_data.C3Cz_GCI)*ones(height(right_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C3C4_GCI)*ones(height(right_data),1)),'Color','blue','LineWidth',2);
plot(1:height(right_data),(mean(right_data.CzC3_GCI)*ones(height(right_data),1)),'Color','green','LineWidth',2);
plot(1:height(right_data),(mean(right_data.CzC4_GCI)*ones(height(right_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4C3_GCI)*ones(height(right_data),1)),'Color',orange,'LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4Cz_GCI)*ones(height(right_data),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(right_data)])
ylim([0 0.2])
title('Cue onset RIGHT')
legend(right_data.Properties.VariableNames{28:33},'Interpreter', 'none')
sgtitle('Granger Causality')
hold off

fig_13=figure(13);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data),left_data.C3Cz_CGCI,'Color','magenta');
hold on
plot(1:height(left_data),left_data.C3C4_CGCI,'Color','blue');
plot(1:height(left_data),left_data.CzC3_CGCI,'Color','green');
plot(1:height(left_data),left_data.CzC4_CGCI,'Color','yellow');
plot(1:height(left_data),left_data.C4C3_CGCI,'Color',orange);
plot(1:height(left_data),left_data.C4Cz_CGCI,'Color','red');
plot(1:height(left_data),(mean(left_data.C3Cz_CGCI)*ones(height(left_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C3C4_CGCI)*ones(height(left_data),1)),'Color','blue','LineWidth',2);
plot(1:height(left_data),(mean(left_data.CzC3_CGCI)*ones(height(left_data),1)),'Color','green','LineWidth',2);
plot(1:height(left_data),(mean(left_data.CzC4_CGCI)*ones(height(left_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4C3_CGCI)*ones(height(left_data),1)),'Color',orange,'LineWidth',2);
plot(1:height(left_data),(mean(left_data.C4Cz_CGCI)*ones(height(left_data),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data)])
ylim([0 0.18])
title('Cue onset LEFT')
legend(left_data.Properties.VariableNames{34:39},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data),right_data.C3Cz_CGCI,'Color','magenta');
hold on
plot(1:height(right_data),right_data.C3C4_CGCI,'Color','blue');
plot(1:height(right_data),right_data.CzC3_CGCI,'Color','green');
plot(1:height(right_data),right_data.CzC4_CGCI,'Color','yellow');
plot(1:height(right_data),right_data.C4C3_CGCI,'Color',orange);
plot(1:height(right_data),right_data.C4Cz_CGCI,'Color','red');
plot(1:height(right_data),(mean(right_data.C3Cz_CGCI)*ones(height(right_data),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C3C4_CGCI)*ones(height(right_data),1)),'Color','blue','LineWidth',2);
plot(1:height(right_data),(mean(right_data.CzC3_CGCI)*ones(height(right_data),1)),'Color','green','LineWidth',2);
plot(1:height(right_data),(mean(right_data.CzC4_CGCI)*ones(height(right_data),1)),'Color','yellow','LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4C3_CGCI)*ones(height(right_data),1)),'Color',orange,'LineWidth',2);
plot(1:height(right_data),(mean(right_data.C4Cz_CGCI)*ones(height(right_data),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(right_data)])
ylim([0 0.18])
title('Cue onset RIGHT')
legend(right_data.Properties.VariableNames{34:39},'Interpreter', 'none')
sgtitle('Conditional Granger Causality')
hold off

% Adding class 0 for left, 1 for right
left_class=array2table(zeros(height(left_data),1),'VariableNames',{'Class'});
right_class=array2table(ones(height(right_data),1),'VariableNames',{'Class'});

left_data=[left_data,left_class];
right_data=[right_data,right_class];

clear left_class right_class

% Merging and shuffling left and right data
dataset=[left_data;right_data];
rng(0); % get the same dataset every time
dataset = dataset(randperm(size(dataset, 1)), :);

% Rescaling data
arr=table2array(dataset(:,1:width(dataset)-1));
colmin=min(arr);
colmax=max(arr);
arr=rescale(arr,'InputMin',colmin,'InputMax',colmax);
dataset(:,1:width(dataset)-1)=array2table(arr);

% Clearing parameters
clear arr colmin colmax C3 Cz C4 check acf lag new_row l_idx r_idx b1 b1_data b2 b2_data cor dif_C3 dif_C4 dif_C5

% Feature Selection
idx1 = fscmrmr(dataset(:,1:width(dataset)-1),dataset.Class);
[idx2,scores] = fscchi2(dataset(:,1:width(dataset)-1),dataset.Class);


% Creating datasets
nf=20;
test_size=50;
dataset1=[dataset(:,idx1(1:nf)),dataset(:,width(dataset))];
dataset2=[dataset(:,idx2(1:nf)),dataset(:,width(dataset))];
% t-test
arr_left=table2array(left_data(:,1:width(left_data)-1));
arr_right=table2array(right_data(:,1:width(right_data)-1));
h=ttest2(arr_left,arr_right);
clear arr_left arr_right
datasett=[dataset(:,h==1),dataset(:,width(dataset))];


% Classification using 39 features
k=14;
q=2;
test_x=dataset(1:test_size,1:width(dataset)-1);
test_y=dataset.Class(1:test_size);
train_x=dataset(test_size+1:height(dataset),1:width(dataset)-1);
train_y=dataset.Class(test_size+1:end);
model=fitcsvm(train_x,train_y,'KernelFunction','polynomial','PolynomialOrder',q,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc_svm=(sum(pred_y==test_y))/test_size;
model=fitctree(train_x,train_y,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc_tree=(sum(pred_y==test_y))/test_size;


% Classification using dataset1
test_x=dataset1(1:test_size,1:nf);
test_y=dataset1.Class(1:test_size);
train_x=dataset1(test_size+1:height(dataset),1:nf);
train_y=dataset1.Class(test_size+1:end);
model=fitcsvm(train_x,train_y,'KernelFunction','polynomial','PolynomialOrder',q,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc1_svm=(sum(pred_y==test_y))/test_size;
model=fitctree(train_x,train_y,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc1_tree=(sum(pred_y==test_y))/test_size;

% Classification using dataset2
test_x=dataset2(1:test_size,1:nf);
test_y=dataset2.Class(1:test_size);
train_x=dataset2(test_size+1:height(dataset),1:nf);
train_y=dataset2.Class(test_size+1:end);
model=fitcsvm(train_x,train_y,'KernelFunction','polynomial','PolynomialOrder',q,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc2_svm=(sum(pred_y==test_y))/test_size;
model=fitctree(train_x,train_y,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc2_tree=(sum(pred_y==test_y))/test_size;

% Classification using datasett
test_x=datasett(1:test_size,1:width(datasett)-1);
test_y=datasett.Class(1:test_size);
train_x=datasett(test_size+1:height(dataset),1:width(datasett)-1);
train_y=datasett.Class(test_size+1:end);
model=fitcsvm(train_x,train_y,'KernelFunction','polynomial','PolynomialOrder',q,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acct_svm=(sum(pred_y==test_y))/test_size;
model=fitctree(train_x,train_y,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acct_tree=(sum(pred_y==test_y))/test_size;






