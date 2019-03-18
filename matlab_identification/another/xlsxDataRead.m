tau1 = xlsread("xls_data/tau_filt1.xlsx");
tau2 = xlsread("xls_data/tau_filt2.xlsx");
tau3 = xlsread("xls_data/tau_filt3.xlsx");

theta1 = xlsread("xls_data/thetaData1.xlsx");
theta2 = xlsread("xls_data/thetaData2.xlsx");
theta3 = xlsread("xls_data/thetaData3.xlsx");

theta_dot1 = xlsread("xls_data/thetaDotData1.xlsx");
theta_dot2 = xlsread("xls_data/thetaDotData2.xlsx");
theta_dot3 = xlsread("xls_data/thetaDotData3.xlsx");

theta_ddot1 = xlsread("xls_data/thetaDDotData1.xlsx");
theta_ddot2 = xlsread("xls_data/thetaDDotData2.xlsx");
theta_ddot3 = xlsread("xls_data/thetaDDotData3.xlsx");

theta = [theta1; theta2; theta3];
theta_dot = [theta_dot1; theta_dot2; theta_dot3];
theta_ddot = [theta_ddot1; theta_ddot2; theta_ddot3];
tau = [tau1; tau2; tau3];
in = [theta, theta_dot, theta_ddot];
out = tau;

in1 = [theta1, theta_dot1, theta_ddot1];

out1 = NNFunc0311(in1);

