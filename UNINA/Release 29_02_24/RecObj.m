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



function f0 = RecObj(x,P_l,P_g,nt,eta_b,eta_i,dt,E_max,E_min,E_b,Pb,Pr_sh,Pr_s,Pr_p,E_0,e_rIII,n_dr_III,n_t_dr3)

P_gt = P_g+x(1:nt)';
P_dr = zeros(nt,1);
for h=1:n_dr_III
    P_dr(1:n_t_dr3,1) = P_dr(1:n_t_dr3,1)+x(nt+(h-1)*n_t_dr3+1:nt+h*n_t_dr3)';
end
P_l = P_l+P_dr;

f0 = -min([P_l';P_gt'])*Pr_sh-P_gt'*Pr_s+P_l'*Pr_p;