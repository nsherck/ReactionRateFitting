function [ Conc, ConcNorm, Mn, MnNorm, Weighted, Pmatrix, MinVal, I, Optimum ] = CGF110215( Mechanism, param, init_c, sample_t, sample_conc, sample_Mn, sample_PDI )
%CGF072915 Utilized in a course grain fitting (CFG) for kinetic data from
%DBU catalyzed ROP of LA or GA. Can be fitted using Mn or Conversion data
%or a weight of either. 
%   Created 8.16.17
%   Revised 8.26.17 to reflect kinetic equations as of 8.26.15 

syms perclb percub L M
cnt = 0 ;
[ri, ci] = size(init_c);
[r, c] = size(param);
[rt,ct] = size(sample_t);
iter = 1; % specifies the resolution
sp = zeros(r, iter); % builds 2D matrix of search values
MW = 144;
for i = 1:r
    sp(i,:) = logspace(log10(param(i,2)), log10(param(i,3)), iter); % buildes the intervals for the param value searches
    %(sp = search space)
end
p = zeros(r,1); % p denotes the parameter and is a vector of parameter values
TotalComb = iter^r
Pmatrix = zeros(TotalComb, r);
i = 1;


%% Defining the Initial Guess

for a = 1:iter 
        p(1) = sp(1,a);
        for b = 1:iter
            p(2) = sp(2,b);
            for c = 1:iter
                p(3) = sp(3,c);
                for d = 1:iter;
                    p(4) = sp(4,d);
                    for e = 1:iter;
                        p(5) = sp(5,e);
                        
                        cnt = cnt + 1 % count
                        
                                          
                        for k = 1:r
                                Pmatrix(cnt, k) = p(k,1);  % creates the matrix for the parameters guessed at the start of each trial
                        end
                        
                        for L =1:ri % iterates through the initial conditions       
                        
                        
                            init = init_c(L,:);
                            tfinal = max(sample_t(:));
                            opts = odeset('RelTol', 1E-7, 'AbsTol', 1E-7);
                            int1 = ode15s(@Mech110215, [0 tfinal], init, opts, p);
                            
                        
                        for j = 1:ct
                            Tf(1,1) = sample_t(L, j);
                            if Tf == 0 
                                ErrConc(j) = 0;
                                ErrMn(j) = 0;
                            else
                                
                                exp_conc = sample_conc(L, j);
                                exp_Mn = sample_Mn(L, j);
                                [m, n] = size(int1.x); % int.x represents the number of steps used by the solver
                        
                                try % checks to see if the function is defined in the region of curiosity
                                    out = deval(int1, Tf);
                                    MnV = MW*( out(16) + out(12) + out(22) + out(26) + out(29) ) / ( out(15) + out(11) + out(21) + out(25) + out(28) );
                                    Mw = MW*(  out(17) + out(13) + out(23) + out(27) + out(30) ) / ( out(16) + out(12) + out(22) + out(26) + out(29) );
                                catch 
                                    out(2) = 0 ; % sets the concentration to 0 
                                    MnV = 0 ;
                                    Mw = 0;
                                end
                               
                                ErrConc(j) = (exp_conc - out(2))^2;
                                ErrMn(j) = ( exp_Mn - (MnV) )^2;
                            end
                        end
                        
                            SSEConc = sum(ErrConc);
                            SSEMn = sum(ErrMn);
                        
                            Conc(cnt,L) = SSEConc;
                            Mn(cnt,L) = SSEMn;
                       
                        for k = 1:r
                            Pmatrix(cnt,k) = p(k,1);
                            
                        end 
                          
                  end  % iterates through initial conditions   
                    end % param 5
                end % param 4
            end % param 3
        end % param 2
end % param 1

%% Concentration
M = 1; 
for M = 1:TotalComb
    Conc(M,ri+1)= sum(Conc(M,1:ri)); % provides the total sum res for conc
end 

L = 1;
for L=1:ri+1
MaxConc = max(Conc(:,L));
ConcNorm(:,L) = Conc(:,L)/MaxConc;
end 

%% Mn
M = 1; 
for M = 1:TotalComb
    Mn(M,ri+1)= sum(Mn(M,2:ri)); % provides the total sum res for conc
end 

L = 1;
for L=1:ri+1
MaxMn = max(Mn(:,L));
MnNorm(:,L) = Mn(:,L)/MaxMn;
end 

%% Weighting

N = 1; 
i = 1;
w1 = 1; % weight for conc.
w2 = 0; % weight for Mn
for N = 1:TotalComb
    Weighted(N,1)= w1*ConcNorm(N,ri+1)+ w2*MnNorm(N,ri+1) ; % provides the total sum res for conc
end 

%% Visualizing Optimum from Weighted Matrix

[MinVal, I] = min(Weighted) ; 
Optimum = Pmatrix(I, :) ; 


for L =1:ri % iterates through the initial conditions    
                        init = init_c(L,:);
                        tfinal = max(sample_t(L,:));
                        opts = odeset('RelTol', 1E-5, 'AbsTol', 1E-6);
                        [T,Y] = ode15s(@Mech110215, [0 tfinal], init, opts, Optimum);
%***********Predicted PDI**************************************%
                        i = 0; 
                        [rp, cp] = size(Y);
                        PDI = zeros(rp, 1);

                        MW = 144; 
                        Mn1 = MW*zeros(rp-1, 1);
                        Mw1 = MW*zeros(rp-1, 1);
                        PDI(1) = 1 ;
                        D = zeros(rp,1);
                        K = zeros(rp,1);
                        C = zeros(rp,1); 
                        
                        D(1) = 0 ;
                        K(1) = 0 ;
                        C = 0 ;

                        for i = 2:rp
    
                            Mn1(i) = MW*( Y(i, 16) + Y(i, 12) + Y(i, 22)+ Y(i, 26) + Y(i,29) ) / ( Y(i, 15) + Y(i,11) + Y(i, 21)+ Y(i, 25) + Y(i,28));
                            Mw1(i) = MW*( Y(i, 17) + Y(i, 13) + Y(i, 23)+ Y(i, 27) + Y(i,30) ) / ( Y(i, 16) + Y(i,12) + Y(i, 22)+ Y(i, 26) + Y(i,29));
                            PDI(i) = Mw1(i)/Mn1(i) ; % Builds Predicted PDI column vector
                            D(i) = Y(i, 21) + Y(i, 25) ; % R terminated chains
                            K(i) = Y(i, 11) + Y(i, 15) ; % Ketene aminal terminated chains
                            C(i) = Y(i, 28) ; % Cyclic Polymer Chains
   
                        end                 
%**************** Visuals *************************************%
                            
                        figure 
                        subplot(2,2,1)
                        set(gcf,'numbertitle','off','name',sprintf('Model - Predicted versus Experimental %2.1f', L))
                        plot(T, Y(:,2), 'r')
                        xlabel('Time [Seconds]')
                        ylabel('[LA]')
                        title('Predicted [LA]')
                        hold on 
                        grid on
                        plot(sample_t(L,:), sample_conc(L,:), 'B*')
                        legend('Predicted [LA]', 'Experimental [LA]')

                        subplot(2,2,2)
                        plot(T, Mn1, 'r')
                        xlabel('Time [Seconds]')
                        ylabel('MW # average')
                        title('Predicted # Average Molecular Weight Profile')
                        hold on 
                        grid on
                        plot(sample_t(L,:), sample_Mn(L,:), 'B*')
                        legend('Predicted Mn', 'Experimental Mn')

                        if sample_PDI(L,1) == 1 
                        [ g, h] = size(sample_PDI(L,:));
                            subplot(2,2,3)
                            plot(T, PDI, 'R')
                            xlabel('Time [Seconds]')
                            ylabel('Poly-dispersity Index')
                            title('Predicted PDI Profile versus Experimental')
                            hold on 
                            grid on
                            plot(sample_t(L,:), sample_PDI(L,2:h), 'B*')
                            legend('Predicted PDI', 'Experimental PDI')%, 'Experimental Values')
                        
                        else
                            subplot(2,2,3)
                            plot(T, PDI, 'R')
                            xlabel('Time [Seconds]')
                            ylabel('Poly-dispersity Index')
                            title('Predicted PDI Profile')
                            grid on
                            legend('Predicted PDI')%, 'Experimental Values')
                        
                        end
                        
                        subplot(2,2,4)
                        plot(T, D(:), 'R')
                        xlabel('Time [Seconds]')
                        ylabel('[Polymer]')
                        title('Profile of [D*+R], [KA*+KA] and [C]')
                        hold on 
                        grid on
                        plot(T, K(:), 'B')
                        plot(T, C(:), 'g')
                        legend('[D*+R]', '[KA*+KA]', '[C]')
    end


end

