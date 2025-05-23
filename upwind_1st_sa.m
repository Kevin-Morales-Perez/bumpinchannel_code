function [Ap_k,Anb_nu_ab_k] = upwind_1st_sa(mf_a,phi_A,Ap_k,Anb_nu_ab_k,mf_t,index)
    %mf_a : mass flux with most influence in face A 
    %phi_A: phi value in adjacent cell 
    %Ap_k: LSE
    %Anb_nu_ab_k: RSE
    %total mass flux across face A

    if mf_a<0 %node A is upwind
    Ap_k(index)=0;
    Anb_nu_ab_k(index)=-mf_t*phi_A;
    else 
    Ap_k(index)=mf_t;
    Anb_nu_ab_k(index)=0;
    end
end
