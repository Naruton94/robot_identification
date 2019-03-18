clc, clearvars, close all
clear dynIdenf

% the number of read data
p0 = 11001;
tsampleo = 0.002;
% read the data from files
fldRead = 'C:\Users\NARUTON\Desktop\20190311\dynIdenf_matlab\20181123_result\DataOutput_Traj1';
[theta, theta_dot, theta_ddot, tau] = loadData(fldRead, p0);
% resample  the data
theta = theta(1:3:p0,:);
tau = tau(1:3:p0,:);
p = size(theta,1);
tsample = 0.006;
% correct the torque of the third axis
tau(:,3) = - tau(:,3);
% scale the torque
Ts = [3.442383e-4*ones(1,3), 1.342773e-4*ones(1,3)];
ratio_trans = 101;
tau = tau * diag(Ts) * ratio_trans;
% the cycle time
T = 20;
% DH parameters given by shanghai dianqi
a = [0, 0, -0.416, -0.4208, 0, 0, 0];
alpha = [0, pi/2, 0, 0, pi/2, -pi/2, 0];
d = [0, 0.1181, 0, 0, 0.1301, 0.1021, 0.0568];
% the gravity
g = [0; 0; -9.80665];
% the priori dynamic parameters
phi_r0 = [2,3.42000000000000,1.26000000000000,0.800000000000000,0.800000000000000,0.350000000000000;
    -0.0165019647020000,0.426439873502640,0.139078355556420,1.39460800000000e-07,-2.63220000000000e-06,8.29500000000000e-11;
    -0.0256289655860000,-0.0124271511264000,-1.07045177400000e-05,0.000225010768000000,-0.000320608410400000,-3.71160097000000e-05;
    -0.0456037864700000,0.386835887459880,0.0411503703017400,-0.00478865628240000,-0.00460317585600000,-0.00683040118320000;
    0.00388830591570456,0.0475012505609916,0.00223362215511204,0.000772782642110361,0.000770846701360872,0.000261196607049744;
    7.63638706167896e-05,3.00901079796083e-06,5.03356734903604e-07,5.04623525623476e-09,-6.16295700192826e-08,-4.06640243700716e-12;
    0.000340020397992622,-0.0528168520144169,-0.00375135880026783,8.09638107508407e-09,-4.41585931650860e-08,1.37724047043447e-11;
    0.00333118131630822,0.131182004246096,0.0264742256951927,0.000537499391466433,0.000531210606854337,0.000261849371147053;
    0.000215211872926742,0.00172637529085978,3.69281505903818e-07,-6.01459794336313e-06,6.03332741029044e-06,-4.02812899656085e-07;
    0.00262696792835856,0.0864328868443469,0.0248796360944246,0.000649992384815521,0.000653231083089124,0.000175309188669306;
    0,0,0,0,0,0;
    0,0,0,0,0,0];
% pars - the input parameters for dynamic identification, in which
% pfilt - number of points to be filtered
% pidenf(1) - number of points used for identification
% pidenf(2) - space of points for identification
% peval(1) - number of points used for torque evaluation
% peval(2) - space of points for torque evaluation
% noise_err - noise level of the elements in matrix K
% cond_max - the maximum condition number allowed for K
% lambda - the weight of the priori dynamic parameters
% fpass - the pass frequency of filter
% tsample - the sample time
% orderfilt - the order of the filter
% pswitch(1) - the multiples of the torque's difference's RMS for points to be
% recoganized as the switch points
% pswitch(2) - the spread of the switching point
% v_zero - the transient velocity of the Coulomb friction
% segErr - the segmentation of the percentage error
pfilt = 3341;
pidenf = [335, 10];
peval = [3341, 1];
noise_err = 1e-6;
cond_max = 100;
lambda = 0;
fpass = 2;
tsample = tsample;
orderfilt = 800;
pswitch = [3, 5, 0];
v_zero = 0.001;
segErr = [0.2 inf];
% parameters of the earlier algorithm of dynamic identification
pars = {a, alpha, d, g, phi_r0, pfilt, pidenf, peval, noise_err, cond_max, lambda, fpass, tsample, orderfilt, v_zero, segErr};
% parameters of the improved algorithm of dynamic identification
% pars = {a, alpha, d, g, phi_r0, pfilt, pidenf, peval, noise_err, cond_max, lambda, fpass, tsample, orderfilt, pswitch, v_zero, segErr};

% number of segments of the percentage error, and the colors in the plots
nSeg = length(segErr);
colors = {'gree','red','red','cyan','magenta','purk'};
% obtain the time of original data and the output dat
t = linspace(0,T/pfilt*p,p)';
teval = linspace(T/pfilt*(p-pfilt)/2,T/pfilt*(p+pfilt)/2,peval(1))';
% expand the error segments by zero
segErr1 = [0.0 segErr];

% ealier identification of dynamic parameters
[phi, tau_pos, tau_pre, tau_filt, errs] = dynIdenf_old(theta, tau, pars);
% improved identification of dynamic parameters
% [phi, tau_pos, tau_pre, tau_filt, errs] = dynIdenf(theta, tau, pars);
xlswrite("tau_pos.xlsx", tau_pos);
% plot the results
figure(1);
for i = 1:6
    subplot(3,2,i);
    hdl = zeros(nSeg+1,1);
    hdl(1) = plot(t,tau(:,i),'color',[0 0 1 0.5],'linewidth',1); hold on
    %plot(teval,tau_filt(:,i),'color',[1 0 0 0.5],'linewidth',2); hold on
    plot(teval,tau_pos(:,i),'-','color','green','linewidth',1,'markersize',0.2); hold on
    for j = nSeg:-1:1
        tau_err = abs((tau_pos(:,i)-tau_filt(:,i)) ./ tau_filt(:,i));
        idx = tau_err > segErr1(j) & tau_err <= segErr1(j+1);
        hdl(j+1) = plot(teval(idx),tau_pos(idx,i),'o','color',[colors{j} 0.5],'linewidth',1,'markersize',1); hold on
    end
    hold off
    xlabel('time / s');
    ylabel('torque / Nm');
    title(['Joint' num2str(i)]);
    leg = legend(hdl,'Experiment',['Predicted - Error<' num2str(segErr*100) '%'],['Predicted - Error>' num2str(segErr*100) '%']);
%     leg = legend('Raw data','Filtered data');
    set(leg,'location','north');
%     xlim([0 24]);
%     ylim([-2 2]);
end