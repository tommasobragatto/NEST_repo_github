%% SOFTWARE FOR PROBABILISTIC OPTIMAL CONTROL OF RESOURCES IN A RENEWABLE ENERGY COMMUNITIY
%  ***********************************************************************************************************************************************
%  RELEASE 29/02/2024 
%  FABIO MOTTOLA AND DANIELA PROTO UNIVERSITY OF NAPLES FEDERICO II, NAPLES (ITALY)
%  fabio.mottola@unina.it, daniela.proto@unina.it
%
%  THE SOFTWARE IS PART OF Project funded under the National Recovery and Resilience Plan (NRRP), 
%  Mission 4 Component 2 Investment 1.3 - Call for tender No. 1561 of 11.10.2022 of Ministero dell Università e della Ricerca (MUR). 
%  European Union – NextGenerationEU. Award Number: Project code PE0000021, Concession Decree No. 1561 of 11.10.2022 
%  adopted by Ministero dell Università e della Ricerca (MUR), CUP E63C22002160007 - Project title “Network 4 Energy Sustainable Transition – NEST.
%  ************************************************************************************************************************************************



function [fi,fu] = RecCon(x,P_l,P_g,nt,eta_b,eta_i,dt,E_max,E_min,E_b,P_b,Pr_sh,Pr_s,Pr_p,E_0,e_rIII,n_dr_III,n_t_dr3)

E_1 =-(x(1:nt)/eta_i.*eta_b.*(x(1:nt)<0)+x(1:nt)/eta_i/eta_b.*(x(1:nt)>0))*dt*P_b/E_b;  % SoC evaluation for all forward time intervals
E_11 = [E_0,E_1];
E = cumsum(E_11);

fi(1:nt,1) = -x(1:nt)'-P_g/eta_i;       % BESS Vs RES generation constraints
fi(nt+1:2*nt,1)=E(2:nt+1)'-E_max;       % BESS maximum capacity
fi(2*nt+1:3*nt,1)=-E(2:nt+1)'+E_min;    % BESS minimum capacity

if n_dr_III>0
    for h=1:n_dr_III
        fi(3*nt+h,1) = -x(nt+(h-1)*n_t_dr3+1:nt+h*n_t_dr3)*ones(n_t_dr3,1)*dt+e_rIII(h,1);       % Energy required by DR Load Type III
    end
end

fu = [];
    