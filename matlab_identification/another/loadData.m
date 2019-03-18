function [theta, theta_dot, theta_ddot, tau] = loadData(fldRead, p)
% map the variables to the file names
filnam_theta = 'DataOutput_theta';
filnam_theta_dot = 'DataOutput_velocity';
filnam_theta_ddot = 'DataOutput_acceleration';
filnam_tau = 'DataOutput_torque';
filnam_phi = 'DataOutput_parameterSet';
filnam_tau_pos = 'DataOutput_torque_pos';
filnam_tau_pre = 'DataOutput_torque_pre';

% set read formats for the variables
datfmtr_theta = 'DataOutput.OutData[%*d]:=%f%*s\n';
datfmtr_theta_dot = 'DataOutput.OutData[%*d]:=%f%*s\n';
datfmtr_theta_ddot = 'DataOutput.OutData[%*d]:=%f%*s\n';
datfmtr_tau = 'DataOutput.OutData[%*d]:=%f%*s\n';
datfmtr_phi = 'DataOutput.OutData[%*d]:=%f%*s\n';
datfmtr_tau_pos = 'DataOutput.OutData[%*d]:=%f%*s\n';
datfmtr_tau_pre = 'DataOutput.OutData[%*d]:=%f%*s\n';

% set write formats for the variables
datfmtw_theta = 'DataOutput.OutData[%d]:=%fF\r\n';
datfmtw_theta_dot = 'DataOutput.OutData[%d]:=%fF\r\n';
datfmtw_theta_ddot = 'DataOutput.OutData[%d]:=%fF\r\n';
datfmtw_tau = 'DataOutput.OutData[%d]:=%fF\r\n';
datfmtw_phi = 'DataOutput.OutData[%d]:=%fF\r\n';
datfmtw_tau_pos = 'DataOutput.OutData[%d]:=%fF\r\n';
datfmtw_tau_pre = 'DataOutput.OutData[%d]:=%fF\r\n';

% set the numbers and dimensions
% ns - the number of sampled points per cycle
% n - the number of joints
n = 6;

% read the data
data = readData(fldRead, {filnam_theta,filnam_theta_dot,filnam_theta_ddot,filnam_tau},...
    [n,n,n,n], {datfmtr_theta,datfmtr_theta_dot,datfmtr_theta_ddot,datfmtr_tau}, [p,p,p,p]);

% output the data
theta = data{1};
theta_dot = data{2};
theta_ddot = data{3};
tau = data{4};
end

function data = readData(fldRead, filnam, ncomp, datfmt, p)
% the function reads data from files
% inputs
% fldRead - the fold path where the files are located, type string
% filnam - the names of files to be read, type cells of strings
% ncomp - the components of variables, type vector
% datfmt - the formats of the data to be read, type cells of strings
% ns - the samples of the data, type vector
% output
% data - the read data

% extension of the files
filext = '.txt';
% the number of variables to be read
nvar = numel(filnam);
% initialize the data output
data = cell(nvar,1);
% open, read and close files
for i = 1:nvar
    data{i} = zeros(p(i), ncomp(i));
    for j = 1:ncomp(i)
        fileID = fopen([fldRead '\' filnam{i} sprintf('%d',j) filext],'r');
        data{i}(:,j) = fscanf(fileID,datfmt{i},p(i));
        fclose(fileID);
    end
end
end