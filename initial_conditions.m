% Loads all initial conditions for orbit determination process

% Initial Epoch
IC.epoch = 2458200.74656; % Julian Date

% Initial State
IC.R_ECI_0 = [6984.45711518852 1612.2547582643 13.0925904314402]'; % km
IC.V_ECI_0 = [-1.67667852227336 7.26143715396544 0.259889857225218]'; % km / s
IC.C_d = 2.0;

% Test Values
% Initial Test Epoch
IC_Test.epoch = 2458150.70833; % Julian Date

% Initial State
IC_Test.R_ECI_0 = [6990.077798814194 1617.465311978378 22.679810569245355]'; % km
IC_Test.V_ECI_0 = [-1.67513972506056 7.27372441330686 0.252688512916741]'; % km / s
IC_Test.C_d = 2.0;