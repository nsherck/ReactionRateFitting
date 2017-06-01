Need to make sure the parameters are set in the code for: (optimized from 9.22.15)

ka1= 4328.8
kd1= 29.15
kp111=ki1= .706
Keq2 = 9.22E-13 (from DFT calculations, Waymouth)
ka2=kd2*Keq2 
kd2= 6.5E16 
kpt= 4.4E10 (proton transfer, right at the limit of being diffusion limited)

For effect of acylation (combination in excess DBU regime):

kac1 = .28


For effect of cyclizaton (self-combination in excess DBU regime):

kac2 = .28


To GET GRAPHS Automatically made, use the CFG110215 fitting function and set iter = 1!!!
Make sure to input dummy variables for the inputs into CFG110215 in order for it to run. 