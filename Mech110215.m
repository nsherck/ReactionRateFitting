function [ dy ] = Mech110215( t, y, param )
%Mech holds kinetic equations in the notation of moments to be numerically integrated.
%   Mech can handle batch and semibatch, both co- and homo-polymerizations.
% The param array needs to be created to correspond to the correct number
% of unknowns to be fitted. NOTE: you cannot have more unknowns then sample
% size. 
%   Created 4.23.2015
%   Revised 6.24.2015 to make the Dact/2. 
%   Revised 7.08.2015 to provide a more comprehensive model, including PDI
%   and a modified termination mechansim to better reflect the revised
%   mechanisms.
%   Revised 7.25.2015 Added Combinatorial Terms for the Re-opening of the
%   Cyclic Polymers. 
%   Revised 7.26.2015 to only work with DBU-Lactide Systems
%   Revised 8.10.2015 to add in exchange and termination reactions
%   Revised 8.17.2015 Termination by kontamination of only DBU, chains made
%   into alcohol and can undergo exhange. 
%   Revised 8.26.2015 Revised the model corresponding to the word document
%   "Kinetic Equations 08.26.15" 
%   Revised 10.08.2015 to include the possibility of deactivation of DBU
%   through the transition to a protonated ketene aminal. 
%   Revised 10.22.2015 to include combination of alcohol with ketene
%   aminal. 
%   Revised 10.26.2015 to include kac2, combination of alcohol terminated
%   KA with itself. Cyclic polymers now present.
%   Revised 10.30.2015 to include termination first order, e.g. activated
%   alcohol goes to alcohol plus DBU.
%   Revised 11.02.2015 to include deactivation of DBU by "acidic"
%   contamination.
%   Finalized 11.12.2015 with derivations performed by hand

dy = zeros(30, 1); % ODEs

%%%%%%%%% Constants with alcohol initiator  %%%%%%%%%%%%%%
Keq1 = 14;  % Alcohol-DBU equilibrium
ka1 = 4320; %774; %param(1); %4328.8; %param(1) ; % initator DBU complexation
kd1 = 29.15; %14/param(1); %16.68;%param(2); %29.15; 
kp111 = .706 ; %2.154;%param(3); %.7055; 




%%%% Constants without alcohol initiator %%%%%%%%%%%%%%%%%%
Keq2 = 9.2202E-13; % DBU-Lactide Equilibrium
ka2 = 6.5E16*Keq2; %2.3357E9*9.2202E-13 ;%param(1)*Keq2;%Keq2*6.5E16;
kd2 = 6.5E16; %2.3357E9 ;%6.5E16; 
kpt = 4.4E10; %6.1585E10 ;% param(2);%4.4E10; %param(4); 
ki1 = .706;%param(2); %2.154; %param(3); %;

%%%% Termination %%%%

kt1 = 0; % termination from acid or contamination (Propagating chain end only, deactivating DBU)
kd = 0; %4.6416; %param(4) ; % Two Active Chains 
kt2 = 0; % Termination of the first order (active chains spontaneously terminate)
kt3 = 0; % Deactivating DBU from acidic contamination (e.g. water)

%%%  Acylation "combination" %%%%%%%%%%%%%

kac1 = .28 ;%param(5); % constant for combination of alcohols with ketene aminal species.
kac2 = 0 ;%param(5); % constant for self acylation KA terminated alcohols



%********************************************************************************************************************%
%********   Rate Equations   ************%

dy(1) = -1*ka1*y(1)*y(3) + kd1*y(4) + -1*ka1*y(1)*( y(25)+y(15) ) + kd1*( y(21)+y(11) ) + -1*ka2*y(1)*y(2) + kd2*y(5) + -1*kpt*y(5)*y(1) + kac1*( (y(11)+y(15))*( y(3)+y(15)+y(25) ) ) + kac2*y(15) + -1*kt3*y(1)*y(7); % DBU 
dy(2) = -1*ki1*y(4)*y(2) + -1*kp111*y(2)*( y(11)+y(21) ) + -1*ka2*y(1)*y(2) + kd2*y(5) + kac2*y(14) ; % monomer
dy(3) = -1*ka1*y(3)*y(1) + kd1*y(4) + -1*kac1*( y(3)*(y(11) + y(15)) ) ; % Initiator (Alcohol Macro-Initiator)
dy(4) = ka1*y(3)*y(1) + -1*kd1*y(4) + -1*ki1*y(2)*y(4) ; % [D-I] DBU MacroInitiator Complex
dy(5) = ka2*y(1)*y(2) + -1*kd2*y(5) + -1*kpt*y(1)*y(5) ; % D1,1 Monomer DBU Complex 
dy(6) = kt1*y(7)*( y(11)+y(21) ) + kd*y(11)*( y(11)+y(21) ) + kd*y(21)*( y(11)+y(21) ) + kt2*( y(11)+y(21) ) + kt3*y(1)*y(7); % Deactivated DBU
dy(7) = -1*kt1*y(7)*( y(11)+y(21) ) + -1*kt3*y(1)*y(7); % Contamination K


%*****************************************************************%
%***************** Excess DBU Pathway ****************************% 


%***************** "Ketene Aminal Active End -DBU-H+ " ***********%
dy(10) = kpt*y(1)*y(5) + -1*kp111*y(2)*y(10) + -1*kt1*y(10)*y(7) + ka1*y(1)*y(14) + -1*kd1*y(10) + -1*kd*y(10)*( y(11)+y(21) ) + -1*kac1*y(10)*( y(3)+y(25)+y(15) ) + -1*kt2*y(10) ; % chain length 1
dy(11) = kpt*y(1)*y(5) + ka1*y(1)*y(15) + -1*kd1*y(11) + -1*kt1*y(11)*y(7) + -1*kd*y(11)*( y(11)+y(21) ) + -1*kac1*y(11)*( y(3)+y(25)+y(15) ) + kac1*y(15)*y(11) + -1*kt2*y(11) ; % zeroeth moment
dy(12) = kpt*y(1)*y(5) + ka1*y(1)*y(16) + -1*kd1*y(12) + -1*kt1*y(12)*y(7) + kp111*y(2)*y(11) + -1*kd*y(12)*( y(11)+y(21) ) + -1*kac1*y(12)*( y(3)+y(25)+y(15) ) + kac1*( y(16)*y(11)+y(15)*y(12) ) + -1*kt2*y(12) ; % first moment
dy(13) = kpt*y(1)*y(5) + ka1*y(1)*y(17) + -1*kd1*y(13) + -1*kt1*y(13)*y(7) + kp111*y(2)*( 2*y(12)+y(11) ) + -1*kd*y(13)*( y(11)+y(21) ) + -1*kac1*y(13)*( y(3)+y(25)+y(15) ) + kac1*( y(17)*y(11)+2*y(16)*y(12)+y(15)*y(13) )  + -1*kt2*y(13); % second moment


%******************** "Ketene Aminal Alcohol Terminated -OH" *********************% 
dy(14) = -ka1*y(1)*y(14) + kd1*y(10) + kt1*y(7)*y(10) + kd*y(10)*( y(11)+y(21) ) + -1*kac1*y(14)*( y(3)+y(25)+2*y(15)+y(11) ) + -1*kac2*y(14)  + kt2*y(10) ; % chain length 1
dy(15) = -ka1*y(1)*y(15) + kd1*y(11) + kt1*y(7)*y(11) + kd*y(11)*( y(11)+y(21) ) + -1*kac1*y(15)*( y(3)+y(25)+2*y(15)+y(11) ) + kac1*y(15)*y(15) + -1*kac2*y(15) + kt2*y(11); % zeroeth moment
dy(16) = -ka1*y(1)*y(16) + kd1*y(12) + kt1*y(7)*y(12) + kd*y(12)*( y(11)+y(21) ) + -1*kac1*y(16)*( y(3)+y(25)+2*y(15)+y(11) ) + 2*kac1*y(16)*y(15) + -1*kac2*y(16)  + kt2*y(12); % first moment
dy(17) = -ka1*y(1)*y(17) + kd1*y(13) + kt1*y(7)*y(13) + kd*y(13)*( y(11)+y(21) ) + -1*kac1*y(17)*( y(3)+y(25)+2*y(15)+y(11) ) + kac1*( 2*y(17)*y(15)+2*y(16)*y(16) ) + -1*kac2*y(17) + kt2*y(13); % second moment


%*****************************************************************%
%********************  Excess Alcohol Pathway ********************%


%*******************  Acitvated Alcohol **************************%

dy(20) = ki1*y(2)*y(4) + ka1*y(1)*y(24) + -1*kd1*y(20) + -1*kp111*y(2)*y(20) + -1*kt1*y(7)*y(20) + -1*kd*y(20)*( y(11)+y(21) ) + kac1*y(10)*y(3)+ -1*kt2*y(20) ; % l = 1 [D-I]
dy(21) = ki1*y(2)*y(4) + ka1*y(1)*y(25) + -1*kd1*y(21) + -1*kt1*y(7)*y(21) + -1*kd*y(21)*( y(11)+y(21) ) + kac1*y(11)*y(3) + kac1*y(25)*y(11) + -1*kt2*y(21)  ; % zeroeth moment
dy(22) = ki1*y(2)*y(4) + ka1*y(1)*y(26) + -1*kd1*y(22) + -1*kt1*y(7)*y(22) + kp111*y(2)*y(21) + -1*kd*y(22)*( y(11)+y(21) ) + kac1*y(12)*y(3) + kac1*( y(26)*y(11)+y(25)*y(12) ) + -1*kt2*y(22); % first moment
dy(23) = ki1*y(2)*y(4) + ka1*y(1)*y(27) + -1*kd1*y(23) + -1*kt1*y(7)*y(23) + kp111*y(2)*( 2*y(22)+y(21) ) + -1*kd*y(23)*( y(11)+y(21) ) + kac1*y(13)*y(3) + kac1*( y(27)*y(11)+2*y(26)*y(12)+y(25)*y(13) ) + -1*kt2*y(23); % second moment

%*******************  Alcohol terminated Polymer Chains *************************%

dy(24) = -1*ka1*y(1)*y(24) + kd1*y(20) + kt1*y(7)*y(20) + kd*y(20)*( y(11)+y(21) ) + -1*kac1*y(24)*( y(11)+y(15) ) + 1*kac1*y(3)*y(14) + kt2*y(20) ; % l = 1 [R1,1]
dy(25) = -1*ka1*y(1)*y(25) + kd1*y(21) + kt1*y(7)*y(21) + kd*y(21)*( y(11)+y(21) ) + -1*kac1*y(25)*( y(11)+y(15) ) + 1*kac1*y(3)*y(15) + kac1*y(25)*y(15) + kt2*y(21) ; % zeroeth moment
dy(26) = -1*ka1*y(1)*y(26) + kd1*y(22) + kt1*y(7)*y(22) + kd*y(22)*( y(11)+y(21) ) + -1*kac1*y(26)*( y(11)+y(15) ) + 1*kac1*y(3)*y(16) + kac1*( y(26)*y(15)+y(25)*y(16) ) + kt2*y(22)  ; % first moment
dy(27) = -1*ka1*y(1)*y(27) + kd1*y(23) + kt1*y(7)*y(23) + kd*y(23)*( y(11)+y(21) ) + -1*kac1*y(27)*( y(11)+y(15) ) + 1*kac1*y(3)*y(17) + kac1*( y(27)*y(15)+2*y(26)*y(16)+y(25)*y(17) ) + kt2*y(23)  ; % second moment

%*****************************************************************%
%********************  Cyclic Polymers ***************************%

dy(28) =  kac2*( y(15) - y(14) ); % zeroeth moment
dy(29) =  kac2*( y(16) - y(14) ); % first moment
dy(30) =  kac2*( y(17) - y(14) ); % second moment



end

