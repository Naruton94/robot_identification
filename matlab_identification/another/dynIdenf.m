function [phi, tau_pos, tau_pre, tau_filt, errs] = dynIdenf(theta, tau, pars)
% the function identify the dynamic parameters

% inputs:
% theta - the raw angular position of joints
% tau - the raw torque of joints

% outputs:
% phi - the identified dynamic parameters
% tau_pos - the torque calculated by currently solved dynamic parameters
% tau_pre - the torque calculated by previously solved dynamic parameters
% errs - the errs of torque, errs = {distErr, filtErr, convErr}
% tau_filt - the measured, filtered and resampled torque of joints
% distErr - the percentage distribution of error in segErr, between tau_pos
% and tau_filt
% filtErr - the error of torque induced by the filtering, between tau_filt and tau
% convErr - convergence error of the iterations, between the tau_pre and
% tau_pos
% rk - the rank of the matrix K, representing the number of dynamic
% parameters that are updated in current iteration

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

% default values
% pfilt = 3341;
% pidenf = [335, 10];
% peval = [3341, 1];
% noise_err = 1e-6;
% cond_max = 100;
% lambda = 0;
% fpass = 2;
% tsample = tsample;
% orderfilt = 800;
% pswitch = [3, 0, 0];
% v_zero = 0.001;
% segErr = [0.2 inf];
% pars = {a, alpha, d, g, phi_r0, pfilt, pidenf, peval, noise_err, cond_max, lambda, fpass, tsample, orderfilt, pswitch, v_zero, segErr};
% phi_r0 = [2,-0.0165019647020000,-0.0256289655860000,-0.0456037864700000,0.00388830591570456,7.63638706167896e-05,0.000340020397992622,0.00333118131630822,0.000215211872926742,0.00262696792835856,0,0;
%     3.42000000000000,0.426439873502640,-0.0124271511264000,0.386835887459880,0.0475012505609916,3.00901079796083e-06,-0.0528168520144169,0.131182004246096,0.00172637529085978,0.0864328868443469,0,0;
%     1.26000000000000,0.139078355556420,-1.07045177400000e-05,0.0411503703017400,0.00223362215511204,5.03356734903604e-07,-0.00375135880026783,0.0264742256951927,3.69281505903818e-07,0.0248796360944246,0,0;
%     0.800000000000000,1.39460800000000e-07,0.000225010768000000,-0.00478865628240000,0.000772782642110361,5.04623525623476e-09,8.09638107508407e-09,0.000537499391466433,-6.01459794336313e-06,0.000649992384815521,0,0;
%     0.800000000000000,-2.63220000000000e-06,-0.000320608410400000,-0.00460317585600000,0.000770846701360872,-6.16295700192826e-08,-4.41585931650860e-08,0.000531210606854337,6.03332741029044e-06,0.000653231083089124,0,0;
%     0.350000000000000,8.29500000000000e-11,-3.71160097000000e-05,-0.00683040118320000,0.000261196607049744,-4.06640243700716e-12,1.37724047043447e-11,0.000261849371147053,-4.02812899656085e-07,0.000175309188669306,0,0]';

% the counter to count the iterations, and store the dynamic parameters in
persistent a alpha d g phi_r0 pfilt pidenf peval noise_err cond_max lambda fpass tsample orderfilt pswitch v_zero segErr
persistent count phi_pre n setResample setEvalTorq m isLevel nparJoint nparMinSet phi_r Rpk Rphi num den
% initialize the parameters
if isempty(count)
    % set the numbers and dimensions
    % n - the number of joints
    % pfilt - the number of remained points after filtering and truncation
    % pidenf - the number of resampled points for parameter identification
    % per cycle, [points, interval]
    % peval - the number of resampled points for torque evaluation per
    % cycle, [points, interval]
    % setResample - the index of sampled points for parameter
    % identification
    % setEvalTorq - the index of sampled points for torque evaluation
    n = size(theta,2);
    pfilt = pars{6};
    pidenf = pars{7};
    peval = pars{8};
    setResample = ((0:pidenf(1)-1)*pidenf(2))+1;
    setEvalTorq = ((0:peval(1)-1)*peval(2))+1;
    
    % calculate the number of dynamic parameters, considering the gravity
    % g - the gravity
    % m - the spatial dimensions
    % isLevel - whether the base is horizontal or declined
    % nparJoint - the number of parameters per joint;
    % nparMinSet - the number of minimum set of parameters, 48(horizontal base)
    % or 50 (declined base)
    g = pars{4};
    m = length(g);
    isLevel = g(1)==0 && g(2)==0;
    nparJoint = (m^2+3*m+6)/2;
    nparMinSet = [nparJoint*n - 3*n - 4, nparJoint*n - 3*n - 6];
    
    % the D-H parameters and Regroup matrix
    a = pars{1};
    alpha = pars{2};
    d = pars{3};
    [Rpk, Rphi] = regroup(isLevel, nparJoint, nparMinSet(isLevel+1), a(1:n+1), alpha(1:n+1), d(1:n+1));
    
    % about the least-square resolution
    % noise_err - the threshold for elements in K not zero
    % cond_max - the maximum condition number of K
    % lambda - the factor for considering the priori
    noise_err = pars{9};
    cond_max = pars{10};
    lambda = pars{11};

    % set parameters of the filter
    % fpass - pass frequency of the filter
    % tsample - the sample time
    % orderfilt - the order of filter
    % num - the numerator of the transfer function of the filter
    % den - the denominator of the transfer function of the filter
    fpass = pars{12};
    tsample = pars{13};
    orderfilt = pars{14};
%     [num, den] = getCoeffs(fpass);
    num = myfir1(orderfilt, fpass*tsample*2);
    den = 1;
    
    % set parameters for detecting switching of the robot joints
    % pswitch - three elements corresponding to the RMS ratio, the adjacent points
    % number at each switch point, the weight for considering these value
    % in identification
    pswitch = pars{15};
    
    % set the parameter for smooth the Coulomb torque when switching occurs
    % v_zero - tau_f = tau_fc * atan(theta_dot / v_zero)
    v_zero = pars{16};
    
    % set the priori values of the parameters
    % phi_r0 the full form of the priori dynamic parameters
    % phi_r - the regrouped form of the priori dynamic parameters
    phi_r0 = pars{5};
    phi_r = Rphi * reshape(phi_r0(1:nparJoint,1:n),n*nparJoint,1);
    
    % set the error segments
    % segErr - segments of the error
    segErr = pars{17};
    
    % the counter
    count = 0;
    % the previous parameters
    phi_pre = phi_r;
end

% differentialize the angular data

theta_dot = diff(theta, 1, 1) / tsample;

theta_ddot = diff(theta_dot, 1, 1) / tsample;

theta_dot = theta_dot(2:end,:);

theta = theta(3:end,:);
tau = tau(3:end,:);

% filter the data and evaluate the filtering error of the torque
[theta, theta_dot, theta_ddot, tau, filtErr] = filtTrimData(theta, theta_dot, theta_ddot, tau, pfilt, num, den);
% theta_dot = diff(theta)/0.006;
% theta_ddot = diff(theta_dot)/0.006;
% theta_dot = [theta_dot(1,:); theta_dot];
% theta_ddot = [theta_ddot(1,:); theta_ddot; theta_ddot(end,:)];
% parameter identification by fordKinematics() and lsqSVD().
% calculate the K according to forward kinematics


% save thetaData theta;
% save thetaDotData theta_dot;
% save thetaDDotData theta_ddot;

xlswrite("thetaData3.xlsx", theta);
xlswrite("thetaDotData3.xlsx", theta_dot);
xlswrite("thetaDDotData3.xlsx", theta_ddot);

K = fordKinematics(theta(setResample,:), theta_dot(setResample,:), theta_ddot(setResample,:), g, a, alpha, d, Rpk, v_zero);

% resolve the dynamic parameters by least-square method with SVD
% decomposition, and record the rank of matrix K
[phi, rk] = lsqSVD(K, tau(setResample,:), phi_pre, count, phi_r, noise_err, cond_max, lambda, pswitch);
% load phi_1
% phi_temp = phi;
% load phi_8
% phi_temp = phi + phi_temp;
% load phi_9
% phi_temp = phi + phi_temp;
% phi = phi_temp / 3;

% evaluate the torque by fordKinematics() and evalTorque()
% calculate the K according to forward kinematics
K = fordKinematics(theta(setEvalTorq,:), theta_dot(setEvalTorq,:), theta_ddot(setEvalTorq,:), g, a, alpha, d, Rpk, v_zero);

% evaluate the torques and evaluate the error of the torque
tau_filt = tau(setEvalTorq,:);
[tau_pos, tau_pre, distErr, convErr] = evalTorque(K, phi, phi_pre, tau_filt, segErr);

% pack and output the errors
errs = {distErr, filtErr, convErr, rk};

% increment counter and update the previous dynamic parameters
count = count + 1;
phi_pre = phi;
end

function [theta, theta_dot, theta_ddot, tau, filtErr] = filtTrimData(theta, theta_dot, theta_ddot, tau,...
    pfilt, num, den)
% this function filter the data

% inputs
% theta - the raw angular position of joints
% theta_dot - the raw angular velocity of joints
% theta_ddot - the raw angular acceleration of joints
% tau - the raw torque of joints

% parameters
% pfilt - the number of remained points after filtering and truncation
% num - numerator of the transfer function of the filter
% den - denominator of the transfer function of the filter

% outputs
% theta - the filtered, resampled angular position of joints
% theta_dot - the filtered, resampled angular velocity of joints
% theta_ddot - the filtered, resampled angular acceleration of joints
% tau- the filtered, resampled torque of joints
% filtErr - the error of torque induced by error

% obtain the dimensions of data
[p, n] = size(theta);
% obtain the index set of the remained points after filtering
setFilt = floor((p-pfilt)/2)+1:floor((p-pfilt)/2)+pfilt;
% initialize the filtering error
filtErr = zeros(1,n);
% filter the data
for i = 1:n
    theta(:,i) = myfiltfilt(num,den,theta(:,i));
    theta_dot(:,i) = myfiltfilt(num,den,theta_dot(:,i));
    theta_ddot(:,i) = myfiltfilt(num,den,theta_ddot(:,i));
    tauf = myfiltfilt(num,den,tau(:,i));
%     tauf = tau(:,i);
    % calculate the filtering error of the torque
    filtErr(i) = rms(tauf(setFilt)-tau(setFilt,i))/rms(tau(setFilt,i));
    tau(:,i) = tauf;
end
% trancate the data
theta = theta(setFilt,:);
theta_dot = theta_dot(setFilt,:);
theta_ddot = theta_ddot(setFilt,:);
tau = tau(setFilt,:);
end

function [Rpk, Rphi] = regroup(isLevel, nparJoint, nparMinSet, a, alpha, d)
% the function calculates the matrix that regroups the dynamic
% parameters into a minimum set
% inputs:

% parameters
% nparJoint - the number dynamic parameters per joint
% nparMinSet - the number of dynamic parameters in minimum set
% a, alpha, d - the D-H parameters

% outputs:
% Rpk - the matrix for regrouping the matrix K
% Rphi - the matrix for regrouping the dynamics parameters phi
% tau = K*phi = (K*Rpk)*(Rphi*phi)

% obtain the number of joints
n = size(a,2)-1;
% initialize Rp1, Rp1 is the complete form of Rpk
Rp1 = eye(nparJoint*n, nparJoint*n);
% the index of parameters to be truncated in each joint or link
% [m, mcx, mcy, mcz, Ioxx, Ioxy, Ioxz, Ioyy, Ioyz, Iozz, cf, cv]
nops = [1 4 8];
% the loop indicates the number of joints
for j = n:-1:2
    % the revised elements in the matrix Rp1
    temp = [-1                                     0                                       0
        -a(j)                                  0                                       0
        d(j+1)*sin(alpha(j))                   sin(alpha(j))                           0
        -d(j+1)*cos(alpha(j))                  -cos(alpha(j))                          0
        -d(j+1)^2                              -2*d(j+1)                               -1
        -a(j)*d(j+1)*sin(alpha(j))             -a(j)*sin(alpha(j))                     0
        a(j)*d(j+1)*cos(alpha(j))              a(j)*cos(alpha(j))                      0
        -(a(j)^2+d(j+1)^2*cos(alpha(j))^2)     -2*d(j+1)*cos(alpha(j))^2               -cos(alpha(j))^2
        -d(j+1)^2*cos(alpha(j))*sin(alpha(j))  -2*d(j+1)*cos(alpha(j))*sin(alpha(j))   -cos(alpha(j))*sin(alpha(j))
        -(a(j)^2+d(j+1)^2*sin(alpha(j))^2)     -2*d(j+1)*sin(alpha(j))^2               -sin(alpha(j))^2
        0                                      0                                       0
        0                                      0                                       0];
    % revise the elements in the matrix Rp1
    Rp1(nparJoint*(j-2)+1:nparJoint*(j-1), nparJoint*(j-1)+nops) = temp;
    Rp1(nparJoint*(j-1)+5, nparJoint*(j-1)+8) = 1;
end
% if the base of the robot is horizontal, the second and third
% dynamic parameters, i.e. mcx, mcy can be truncated
if isLevel ~= 0
    fnops = [1 2 3 4 5 6 7 8 9];
    % otherwise, mcx, mcy must be considered
else
    fnops = [1 4 5 6 7 8 9];
end
% the number of truncated parameters in the first joint or link
nop = numel(nops);
% the number of truncated parameters in each of other joint
fnop = numel(fnops);
% initialize the index set of all truncated parameters
nopset = zeros(1, fnop+nop*(n-1));
% assign the index set of all truncated parameters
nopset(1:fnop) = fnops;
for j = 2:n
    nopset(fnop+nop*(j-2)+1:fnop+nop*(j-1)) = nparJoint*(j-1)+nops;
end
% the index of the minimum set of dynamic parameters
pset = setdiff(1:nparJoint*n, nopset);
% extract the pset columns of Rp1 and assign them to Rpk
Rpk = zeros(nparJoint*n,nparMinSet);
Rpk(:,:) = Rp1(:,pset);
% inverse of Rp1
invRp1 = eye(nparJoint*n, nparJoint*n) / Rp1;
% extract the pset rows of invRp1 and assign them to Rphi
Rphi = invRp1(pset,:);
end

function K = fordKinematics(theta, theta_dot, theta_ddot,...
    g, a, alpha, d, Rpk, v_zero)
% this function calculates the K kinematics, noting that the
% modified D-H model is used

% inputs:
% theta - the angular position of joints
% theta_dot - the angular velocities of joints
% theta_ddot - the angular acceleration of joints

% parameters:
% g - the gravity
% a, alpha, d - the D-H parameters
% Rpk - the matrix for regrouping K
% v_zero - the transient velocity for smooth the Coulomb friction torque

% outputs
% K - the matrix that makes tau = K*phi.

% obtain dimension of data
[p, n] = size(theta);
m = length(g);
nparJoint = (m^2+3*m+6)/2;
% initialize K
K = zeros(n*p,n*nparJoint);
% calculate the K
% the outer-most loop indicates the number of sample points
for k = 1:p
    % the transformation matrix from reference i to reference 0,
    % and its reverse
    o_X_i = zeros(2*m,2*m,n);
    i_X_o = zeros(2*m,2*m);
    % the absolute velocity and acceleration of reference i with components in
    % itself
    i_V_i = zeros(2*m,1);
    i_A_i = zeros(2*m,1);
    % the secon loop indicates the number of joints
    for i = 1:n
        % the rotation matrix from reference i to i-1
        i1_R_i = [cos(theta(k,i)), -sin(theta(k,i)), 0;...
            sin(theta(k,i))*cos(alpha(i)), cos(theta(k,i))*cos(alpha(i)), -sin(alpha(i));...
            sin(theta(k,i))*sin(alpha(i)), cos(theta(k,i))*sin(alpha(i)), cos(alpha(i))];
        % the offset vector from reference i to i-1
        i1_p_i = [a(i); -d(i+1)*sin(alpha(i)); d(i+1)*cos(alpha(i))];
        % the spatial form of the homogeneous transformation
        % matrix, and its reverse
        i1_X_i = [i1_R_i, zeros(m,m); skew(i1_p_i)*i1_R_i, i1_R_i];
        i_X_i1 = [i1_R_i', zeros(m,m); i1_R_i'*skew(i1_p_i)', i1_R_i'];
        % the relative velocity and acceleration of reference
        % i to i-1 with components in reference i
        i_V_i1_i = [zeros(m-1,1); theta_dot(k,i); zeros(m,1)];
        i_A_i1_i = [zeros(m-1,1); theta_ddot(k,i); zeros(m,1)];
        % if it is the frist joint, the abolute transformation,
        % velocity and acceleration are equal to the relative ones
        if i == 1
            o_X_i(:,:,i) = i1_X_i;
            i_X_o = i_X_i1;
            i_V_i = i_V_i1_i;
            i_A_i = i_A_i1_i;
            % otherwise, obtain the absolute transformation,
            % velocity and acceleration by iteration
        else
            o_X_i(:,:,i) = o_X_i(:,:,i-1) * i1_X_i;
            i_X_o = i_X_i1 * i_X_o;
            i_V_i1 = i_X_i1 * i_V_i;
            i_V_i = i_V_i1 + i_V_i1_i;
            i_A_i = i_X_i1 * i_A_i - skew6(i_V_i1_i) * i_V_i1 + i_A_i1_i;
        end
        % the absolute angular velocity and acceleration (different from classical one) of reference i with
        % components in itself
        i_omega_i = i_V_i(1:m);
        i_omega_dot_i = i_A_i(1:m);
        % the classical absolute acceleration of reference i with
        % components in itself
        i_d_ddot_i = i_A_i(m+1:2*m) + skew(i_V_i(1:m)) * i_V_i(m+1:2*m) - i_X_o(1:m,1:m) * g;
        % i_Am_i is a intermediate matrix, please refer to the book
        % Springer Hand of Robotics,
        i_Am_i = [zeros(m,1),-skew(i_d_ddot_i),lin(i_omega_dot_i)+skew(i_omega_i)*lin(i_omega_i);...
            i_d_ddot_i,skew(i_omega_dot_i)+skew(i_omega_i)^2,zeros(m,2*m)];
        % calculate the matrix K by iterations
        for j = 1:i
            temp = (i_X_o * o_X_i(:,:,j))' * i_Am_i;
            K(n*(k-1)+j,nparJoint*(i-1)+1:nparJoint*i) = [temp(m,:), (i==j)*atan(theta_dot(k,i)/v_zero), (i==j)*theta_dot(k,i)];
        end
    end
end
% regroup the K corresponding to the minimum set of dynamic
% parameters
K = K * Rpk;

% transform a R3 vector into a skew matrix
    function y = skew(x)
        y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0];
    end
% transform a R6 vector into a skew matrix
    function y = skew6(x)
        y = [skew(x(1:3)), zeros(3,3); skew(x(4:6)), skew(x(1:3))];
    end
% transform the angular velocity and acceleration vector into a
% matrix for constructing the linear form of dynamics
    function y = lin(x)
        y = [x(1),x(2),x(3),0,0,0;...
            0,x(1),0,x(2),x(3),0;...
            0,0,x(1),0,x(2),x(3)];
    end
end

function [phi, rk] = lsqSVD(K, tau, phi_pre, count,...
    phi_r, noise_err, cond_max, lambda, pswitch)
% the function solves the least-square with SVD decomposition
% inputs:
% K - the matrix calculated from the measured motion states, the D-H
% parameters and the gravity
% tau - the measured torque
% phi_pre - the previously solved dynamic parameters
% count - the counter

% parameters:
% phi_r - the priori values of the dynamic parameters
% noise_err - the threshold for elements in K not zero
% cond_max - the maximum condition number of K
% lambda - the factor for considering the priori
% pswitch - the parameters for considering the switch points

% outputs:
% phi - the solved dynamic parameters in current iteration
% rk - the rank of the matrix K

% obtain the dimensions of data
[pidenf, n] = size(tau);
nparMinSet = length(phi_pre);
% calculate the step_ratio according to the counter
step_ratio = 1 / (count + 1);
% calculate the RMS of tau, and reshape it into a vector for scaling
tau_rms = rms(tau);
scale_tau = repmat(tau_rms', size(tau,1), 1);
% scale tau and K by the RMS factor
tau = bsxfun(@rdivide, tau, tau_rms);
K = bsxfun(@rdivide, K, scale_tau);
% calculate the RMS of diff tau, and rehape it into a vector for scaling
tau_rms = bsxfun(@rdivide, abs(diff(tau)), rms(diff(tau)));
tau_rms = [tau_rms(1,:); tau_rms];
tau_rms = tau_rms > pswitch(1);
[row, col] = find(tau_rms);
for i = 1:numel(row)
    tau_rms(min(max(row(i)-pswitch(2):row(i)+pswitch(2),1),size(tau_rms,1)), col(i)) = true;
end
tau_rms = tau_rms' * pswitch(3) + (~tau_rms') * 1;
scale_tau = tau_rms(:);
% scale tau and K again
tau = bsxfun(@times, tau, tau_rms');
K = bsxfun(@times, K, scale_tau);
% reshape tau into a vector
tau = tau';
tau = tau(:);
% the scale factor of tau by the absolute value
% scale_tau = abs(tau).^0;
% scale tau and K by the absolute value factor
% tau = tau ./ scale_tau;
% K = K ./ scale_tau;
% detect near-zero elements in K and set them to be zero
K = K .* (abs(K)>noise_err);
% calculate the norm of each column of K
scv = ((vecnorm(K)==0) + vecnorm(K))';
% normalize each column of K and combine it with a lambda-scaled
% identity
K = [K ./ repmat(scv',n*pidenf,1); lambda * eye(nparMinSet)];
% Initialize the current dynamic parameters by the previous dynamic parameters
% and scale it by the norm vector of K
phi = phi_pre .* scv;
% combine the torque matrix with the priori dynamic parameters
tau = [tau; lambda * phi_r.*scv];
% SVD decomposition of K, i.e. K = U*S*V'
[U,S,V] = svd(K);
% obtain the singular values
s = diag(S);
% truncate the smaller singular values to make the condition number
% below cond_max
e1 = s >= max(s)/cond_max;
% calculate the rank of the truncated matrix K
rk = length(find(e1));
% transform the torque by U
tau = U' * tau;
% transform the dynamic parameters by V
phi = V' * phi;
% get the solution in the transoformed form
phi(e1) = tau(e1) ./ s(e1);
% restore the solved dynamic parameters by inverse transformation
% and scaling
phi = V * phi;
phi = phi ./ scv;
% weight the currently calculated dynamic parameters and the
% previous one
phi = phi_pre + step_ratio * (phi - phi_pre);
end

function y_out = myfiltfilt(b,a,x)
% the function filter the data x with the filter specified by coefficients
% b and a

% inputs:
% x - the raw data

% parameters
% b - the numerator of the transfer function of the filter
% a - the denominator of the transfer function of the filter

% outputs
% y - the filtered data

% length of the data
len = size(x,1);
% convert the coefficients into row vector
b = b(:).';
a = a(:).';
% number of elements in the numerator and denominator
nb = length(b);
na = length(a);
% the maximum number of elements
nfilt = max(nb,na);
% length of edge transients
nfact = 3*(nfilt-1);
% set up filter's initial conditions to remove dc offset problems at the
% beginning and end of the sequence
if nb < nfilt
    b = [b, zeros(1,nfilt-nb)];
end
% zero-pad if necessary
if na < nfilt
    a1 = [a, zeros(1,nfilt-na)];
else
    a1 = a;
end
% use sparse matrix to solve system of linear equations for initial conditions
rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
data = [1+a1(2) a1(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
sp = sparse(rows,cols,data);
% zi are the steady-state states of the filter b(z)/a(z) in the state-space
% non-sparse:
% zi = ( eye(nfilt-1) - [-a(2:nfilt).' [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
%      ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
zi = full(sp) \ ( b(2:nfilt).' - a1(2:nfilt).'*b(1) );
% Extrapolate beginning and end of data sequence using a "reflection
% method".  Slopes of original and extrapolated sequences match at
% the end points. This reduces end effects.
y = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)];
% filter, reverse data, filter again, and reverse data again
y = filter(b,a1,y,zi*y(1));
y = y(length(y):-1:1);
y = filter(b,a1,y,zi*y(1));
y = y(length(y):-1:1);
% remove extrapolated pieces of y
y_out = y([nfact+1:len+nfact len+2*nfact+1:length(y)]);
end

function [tau_pos, tau_pre, distErr, convErr] = evalTorque(K, phi, phi_pre, tau_filt, segErr)
% this function evaluate the torques predicted by the model and evaluate
% the error of torque

% inputs:
% K - the K matrix
% phi - the currently resolved dynamics paramters
% phi_pre - the previously resolved dynamics paramters
% tau_filt - the filtered torque
% segErr - the segmentation of the percentage error

% parameters:
% n - the number of joints

% outputs:
% tau_pre - the torque calculated by previously solved dynamic parameters
% tau_pos - the torque calculated by currently solved dynamic parameters
% distErr - the distribution of error
% convErr - the convergence error

% obtain the number of joints and the number of error segmentations
[peval, n] = size(tau_filt);
nSeg = length(segErr);
% calculate tau_pos
tau_pos = K * phi;
tau_pos = reshape(tau_pos,n,peval)';
% calculate tau_pre
tau_pre = K * phi_pre;
tau_pre = reshape(tau_pre,n,peval)';
% calculate the distribution of prediction error
distErr = zeros(nSeg,n);
tau_err = sort(abs((tau_pos - tau_filt)./tau_filt),'descend');
for i = 1:nSeg
    tau_flag = tau_err <= segErr(i);
    for j = 1:n
        idx = sum(find(tau_flag(:,j),1,'first'));
        if isempty(idx) == 0
            distErr(i,j) = (peval - idx + (idx~=peval)) / peval;
        end
    end
end
% calculate the convergence error
tau_err_pos = rms(tau_pos - tau_filt) ./ rms(tau_filt);
tau_err_pre = rms(tau_pre - tau_filt) ./ rms(tau_filt);
convErr = tau_err_pre - tau_err_pos;
end

function bb = myfir1(N,Wn)
% this function calculate the coefficients of the filter

% inputs
% N - the order of the filter
% Wn - the normalized stop frequency

% outputs
% bb - the coefficients of the filter

N = N+1;
odd = mod(N,2);
wind = 0.54-0.46*cos(2*pi*(0:N-1)'/(N-1));
f1 = Wn / 2.0;
c1 = f1;
nhlf = ((N+1) - mod(N+1,2)) / 2;
i1 = odd + 1;
b = zeros(1,nhlf);
if odd
    b(1) = 2 * c1;
end
xn = (0:nhlf-1) + 0.5 * (1 - odd);
c = pi * xn;
c3 = 2 * c1 * c;
b(i1:nhlf)=(sin(c3(i1:nhlf))./c(i1:nhlf));
bb = real([b(nhlf:-1:i1) b(1:nhlf)].*wind(:)');
gain = abs(sum(bb));
bb = bb / gain;
end