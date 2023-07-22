left_data1=left_data(:,1:27);
tab=left_data.C3Cz_GCI./left_data.CzC3_GCI;
tab=array2table(tab,'VariableNames',{'C3CzdivCzC3_GCI'});
left_data1=[left_data1,tab];
tab=left_data.C3Cz_CGCI./left_data.CzC3_CGCI;
tab=array2table(tab,'VariableNames',{'C3CzdivCzC3_CGCI'});
left_data1=[left_data1,tab];
tab=left_data.C3C4_GCI./left_data.C4C3_GCI;
tab=array2table(tab,'VariableNames',{'C3C4divC4C3_GCI'});
left_data1=[left_data1,tab];
tab=left_data.C3C4_CGCI./left_data.C4C3_CGCI;
tab=array2table(tab,'VariableNames',{'C3C4divC4C3_CGCI'});
left_data1=[left_data1,tab];
tab=left_data.CzC4_GCI./left_data.C4Cz_GCI;
tab=array2table(tab,'VariableNames',{'CzC4divC4Cz_GCI'});
left_data1=[left_data1,tab];
tab=left_data.CzC4_CGCI./left_data.C4Cz_CGCI;
tab=array2table(tab,'VariableNames',{'CzC4divC4Cz_CGCI'});
left_data1=[left_data1,tab,left_data(:,40:width(left_data))];

right_data1=right_data(:,1:27);
tab=right_data.C3Cz_GCI./right_data.CzC3_GCI;
tab=array2table(tab,'VariableNames',{'C3CzdivCzC3_GCI'});
right_data1=[right_data1,tab];
tab=right_data.C3Cz_CGCI./right_data.CzC3_CGCI;
tab=array2table(tab,'VariableNames',{'C3CzdivCzC3_CGCI'});
right_data1=[right_data1,tab];
tab=right_data.C3C4_GCI./right_data.C4C3_GCI;
tab=array2table(tab,'VariableNames',{'C3C4divC4C3_GCI'});
right_data1=[right_data1,tab];
tab=right_data.C3C4_CGCI./right_data.C4C3_CGCI;
tab=array2table(tab,'VariableNames',{'C3C4divC4C3_CGCI'});
right_data1=[right_data1,tab];
tab=right_data.CzC4_GCI./right_data.C4Cz_GCI;
tab=array2table(tab,'VariableNames',{'CzC4divC4Cz_GCI'});
right_data1=[right_data1,tab];
tab=right_data.CzC4_CGCI./right_data.C4Cz_CGCI;
tab=array2table(tab,'VariableNames',{'CzC4divC4Cz_CGCI'});
right_data1=[right_data1,tab,right_data(:,40:width(right_data))];

orange=[1,0.7,0];
fig_14=figure(14);
t = tiledlayout(1,2,'TileSpacing','Compact');
subplot(1,2,1)
plot(1:height(left_data1),left_data1.C3CzdivCzC3_GCI,'Color','magenta');
hold on
plot(1:height(left_data1),left_data1.C3CzdivCzC3_CGCI,'Color','blue');
plot(1:height(left_data1),left_data1.C3C4divC4C3_GCI,'Color','green');
plot(1:height(left_data1),left_data1.C3C4divC4C3_CGCI,'Color','yellow');
plot(1:height(left_data1),left_data1.CzC4divC4Cz_GCI,'Color',orange);
plot(1:height(left_data1),left_data1.CzC4divC4Cz_CGCI,'Color','red');
plot(1:height(left_data1),(mean(left_data1.C3CzdivCzC3_GCI)*ones(height(left_data1),1)),'Color','magenta','LineWidth',2);
plot(1:height(left_data1),(mean(left_data1.C3CzdivCzC3_CGCI)*ones(height(left_data1),1)),'Color','blue','LineWidth',2);
plot(1:height(left_data1),(mean(left_data1.C3C4divC4C3_GCI)*ones(height(left_data1),1)),'Color','green','LineWidth',2);
plot(1:height(left_data1),(mean(left_data1.C3C4divC4C3_CGCI)*ones(height(left_data1),1)),'Color','yellow','LineWidth',2);
plot(1:height(left_data1),(mean(left_data1.CzC4divC4Cz_GCI)*ones(height(left_data1),1)),'Color',orange,'LineWidth',2);
plot(1:height(left_data1),(mean(left_data1.CzC4divC4Cz_CGCI)*ones(height(left_data1),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(left_data1)])
ylim([0 12])
title('Cue onset LEFT')
legend(left_data1.Properties.VariableNames{28:33},'Interpreter', 'none')
hold off
subplot(1,2,2)
plot(1:height(right_data1),right_data1.C3CzdivCzC3_GCI,'Color','magenta');
hold on
plot(1:height(right_data1),right_data1.C3CzdivCzC3_CGCI,'Color','blue');
plot(1:height(right_data1),right_data1.C3C4divC4C3_GCI,'Color','green');
plot(1:height(right_data1),right_data1.C3C4divC4C3_CGCI,'Color','yellow');
plot(1:height(right_data1),right_data1.CzC4divC4Cz_GCI,'Color',orange);
plot(1:height(right_data1),right_data1.CzC4divC4Cz_CGCI,'Color','red');
plot(1:height(right_data1),(mean(right_data1.C3CzdivCzC3_GCI)*ones(height(right_data1),1)),'Color','magenta','LineWidth',2);
plot(1:height(right_data1),(mean(right_data1.C3CzdivCzC3_CGCI)*ones(height(right_data1),1)),'Color','blue','LineWidth',2);
plot(1:height(right_data1),(mean(right_data1.C3C4divC4C3_GCI)*ones(height(right_data1),1)),'Color','green','LineWidth',2);
plot(1:height(right_data1),(mean(right_data1.C3C4divC4C3_CGCI)*ones(height(right_data1),1)),'Color','yellow','LineWidth',2);
plot(1:height(right_data1),(mean(right_data1.CzC4divC4Cz_GCI)*ones(height(right_data1),1)),'Color',orange,'LineWidth',2);
plot(1:height(right_data1),(mean(right_data1.CzC4divC4Cz_CGCI)*ones(height(right_data1),1)),'Color','red','LineWidth',2);
xlabel('Sample index')
xlim([1 height(right_data1)])
ylim([0 12])
title('Cue onset RIGHT')
legend(right_data1.Properties.VariableNames{28:33},'Interpreter', 'none')
sgtitle('Granger Causality and Conditional Granger Causality')
hold off

dataset3=[left_data1;right_data1];
rng(0); % get the same dataset every time
dataset3 = dataset3(randperm(size(dataset3, 1)), :);

arr=table2array(dataset3(:,1:33));
colmin=min(arr);
colmax=max(arr);
arr=rescale(arr,'InputMin',colmin,'InputMax',colmax);
dataset3(:,1:33)=array2table(arr);

% Feature Selection
idx1 = fscmrmr(dataset3(:,1:width(dataset3)-1),dataset3.Class);
idx2 = fscchi2(dataset3(:,1:width(dataset3)-1),dataset3.Class);

% t-test
arr_left=table2array(left_data1(:,1:width(left_data1)-1));
arr_right=table2array(right_data1(:,1:width(right_data1)-1));
h=ttest2(arr_left,arr_right);
clear arr_left arr_right
datasett1=[dataset3(:,h==1),dataset3(:,width(dataset3))];


% Creating 2 datasets
nf=20;
test_size=50;
dataset4=[dataset3(:,idx1(1:nf)),dataset3(:,width(dataset3))];
dataset5=[dataset3(:,idx2(1:nf)),dataset3(:,width(dataset3))];

% Classification using 33 features
k=14;
q=2;
test_x=dataset3(1:test_size,1:width(dataset3)-1);
test_y=dataset3.Class(1:test_size);
train_x=dataset3(test_size+1:height(dataset3),1:width(dataset3)-1);
train_y=dataset3.Class(test_size+1:end);
model=fitcsvm(train_x,train_y,'KernelFunction','polynomial','PolynomialOrder',q,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc3_svm=(sum(pred_y==test_y))/test_size;
model=fitctree(train_x,train_y,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc3_tree=(sum(pred_y==test_y))/test_size;


% Classification using dataset4
test_x=dataset4(1:test_size,1:nf);
test_y=dataset4.Class(1:test_size);
train_x=dataset4(test_size+1:height(dataset4),1:nf);
train_y=dataset4.Class(test_size+1:end);
model=fitcsvm(train_x,train_y,'KernelFunction','polynomial','PolynomialOrder',q,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc4_svm=(sum(pred_y==test_y))/test_size;
model=fitctree(train_x,train_y,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc4_tree=(sum(pred_y==test_y))/test_size;

% Classification using dataset5
test_x=dataset5(1:test_size,1:nf);
test_y=dataset5.Class(1:test_size);
train_x=dataset5(test_size+1:height(dataset5),1:nf);
train_y=dataset5.Class(test_size+1:end);
model=fitcsvm(train_x,train_y,'KernelFunction','polynomial','PolynomialOrder',q,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc5_svm=(sum(pred_y==test_y))/test_size;
model=fitctree(train_x,train_y,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acc5_tree=(sum(pred_y==test_y))/test_size;

% Classification using datasett1
test_x=datasett1(1:test_size,1:width(datasett1)-1);
test_y=datasett1.Class(1:test_size);
train_x=datasett1(test_size+1:height(datasett1),1:width(datasett1)-1);
train_y=datasett1.Class(test_size+1:end);
model=fitcsvm(train_x,train_y,'KernelFunction','polynomial','PolynomialOrder',q,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acct1_svm=(sum(pred_y==test_y))/test_size;
model=fitctree(train_x,train_y,'CrossVal','on','KFold',k);
model=model.Trained{1};
pred_y=predict(model,test_x);
acct1_tree=(sum(pred_y==test_y))/test_size;
