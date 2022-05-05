import numpy as np
import numpy.matlib
    
def CTM(CTM_param, rho, supp_p, D_on_ramp_s_str, s_s_str, r_s_str, r_s_ctrl_str, l_str, e_str, intersection_dim, E_cal_out, phi_1, k, opt): 
# + CTM_param: structure with the number of parameters
#    fields->
#        - len (km)[Nx1] : the length of the N cells
#        - v_bar (km/h)[Nx1] : cells free-flow speed
#        - w (km/h)[Nx1] : congestion wave speed
#        - q_max (veh/h)[Nx1] : maximum cell capacity
#        - rho_max (veh/km) [Nx1] : maximum cell density
#        - T (h) [Nx1]: length of each interval
#        - N : number of cells
#        - supply_N_plus [1xK]: the supply of cell N+1 (it is a fictitious value that is used to make the model fit the data properly)
# + rho [Nx1] (veh/km): cell density during the interval k
# + phi_1 (veh/h) : the in-flow of the first cell
# + r_s_str (veh/h) : road-to-station flow
# + s_s_str (veh/h) : station-to-road flow
# + k : time iteration
    
# Output:
#        - rho_plus [N+1x1] (veh/km): the density at the time interval k+1
#        - phi [Nx1] (veh/h): the flow between the cells
    
    try:
        ## Parameters
        N = CTM_param.N
        beta_s_str = CTM_param.beta_s_str
        r_s_max_str = CTM_param.r_s_max_str
        beta = CTM_param.beta
        D_on_ramp = CTM_param.D_on_ramp
        phi_priority = CTM_param.phi_priority
        r_max = CTM_param.r_max
        v_bar = CTM_param.v_bar
        w = CTM_param.w
        q_max = CTM_param.q_max
        rho_max = CTM_param.rho_max
        len_ = CTM_param.len
        T = CTM_param.T
        supply_N_plus = CTM_param.supply_N_plus(k)
        
        ## Field names of structures
        field_beta_s = string(fieldnames(beta_s_str))
        field_r_s_ctrl = string(fieldnames(r_s_ctrl_str))
        field_s_s = string(fieldnames(s_s_str))
        field_e = string(fieldnames(e_str))
        field_r_s_max = string(fieldnames(r_s_max_str))
        field_D_on_ramp_s = string(fieldnames(D_on_ramp_s_str))
        field_r_s = string(fieldnames(r_s_str))
        field_l = string(fieldnames(l_str))
        field_E_cal_out = string(fieldnames(E_cal_out))
        
        ## Demand and Supply of cells
        D = np.zeros((N,1))
        S = np.zeros((N,1))
        for ell in np.arange(1,N+1).reshape(-1):
            # demand ell
            tmp_demand_ell = np.array([(1 - beta(ell) - sum(getattr(beta_s_str,(field_beta_s(ell))))) * v_bar(ell) * rho(ell),q_max(ell)])
            D[ell] = np.amin(tmp_demand_ell)
            # supply ell
            tmp_supply_ell = np.array([w(ell) * (rho_max(ell) - rho(ell)),q_max(ell)])
            S[ell] = np.amin(tmp_supply_ell)
        
        ## Demand of service stations on ramps
        for ell in np.arange(1,N+1).reshape(-1):
            for i in np.arange(1,supp_p(ell,3)+1).reshape(-1):
                tmp_D_on_ramp_s = np.array([getattr(r_s_ctrl_str,(field_r_s_ctrl(ell)))(i,k) + getattr(e_str,(field_e(ell)))(i),getattr(r_s_max_str,(field_r_s_max(ell)))(i)])
                getattr[D_on_ramp_s_str,[field_D_on_ramp_s[ell]]][i] = np.amin(tmp_D_on_ramp_s)
        
        ## Flow entering cell ell from ell-1
        phi,r_s_str = compute_phi_on_ramp(D,D_on_ramp_s_str,r_s_str,S,phi_priority,E_cal_out,phi_1(k),supp_p,N)
        phi[N + 1] = np.amin(np.array([D(N),supply_N_plus]))
        
        ## In and Out flow
        
        # Phi_minus(ell) : outflow of cell ell
        Phi_minus = np.zeros((N,1))
        
        # Phi_plus(ell) : inflow of cell ell
        Phi_plus = np.zeros((N,1))
        
        # r(ell) : flow of on-ramp ell
        r = np.zeros((N,1))
        
        # s(ell) : flow of off-ramp ell
        s = np.zeros((N,1))
        
        # s_s(ell) : flow of off-ramp service station (ell,.)
        s_s = np.zeros((N,1))
        sum_intersection_dim = np.sum(intersection_dim, 1-1)
        
        ## Computation of the components of s_s_str
        for ell in np.arange(1,N+1).reshape(-1):
            for i in np.arange(1,supp_p(ell,2)+1).reshape(-1):
                getattr[s_s_str,[field_s_s[ell]]][i] = getattr(beta_s_str,(field_beta_s(ell)))(i) / (1 - getattr(beta_s_str,(field_beta_s(ell)))(i)) * phi(ell + 1)
        
        # Computation of Phi_plus, s_s_str, s and Phi_minus
        for ell in np.arange(1,N+1).reshape(-1):
            Phi_plus[ell] = phi(ell) + r(ell) + sum(getattr(r_s_str,(field_r_s(ell))))
            s_s[ell] = sum(getattr(s_s_str,(field_s_s(ell))))
            s[ell] = phi(ell + 1) * (beta(ell) / (1 - beta(ell)))
            Phi_minus[ell] = phi(ell + 1) + s(ell) + s_s(ell)
        
        # Sanity check (flows have to be positive)
        if not np.all(Phi_minus >= 0) :
            fmt = np.array(['The vector Phi_minus is: [',np.matlib.repmat('%g, ',1,np.asarray(Phi_minus).size - 1),'%g]\n'])
            fmt.write(Phi_minus % ())
            keyboard
        
        ######################################
        
        if not np.all(Phi_plus >= 0) :
            fmt = np.array(['The vector Phi_plus is: [',np.matlib.repmat('%g, ',1,np.asarray(Phi_plus).size - 1),'%g]\n'])
            fmt.write(Phi_plus % ())
        
        ## Vehicles waiting to leave the charging station
        for ell in np.arange(1,N+1).reshape(-1):
            if sum_intersection_dim(ell) != 0:
                setattr(e_plus_str,field_e(ell),np.zeros((sum_intersection_dim(ell),1)))
                setattr(e_plus_str,field_e(ell),getattr(e_str,(field_e(ell))) + getattr(r_s_ctrl_str,(field_r_s_ctrl(ell)))(:,k) - getattr(r_s_str,(field_r_s(ell))))
            else:
                setattr(e_plus_str,field_e(ell),[])
        
        ## Queue lenght
        
        # l_plus_str (veh):  Queue lenght in the interval k+1
        c2 = np.ones((N,1))
        c3 = np.ones((N,1))
        for ell in np.arange(1,N+1).reshape(-1):
            if sum_intersection_dim(ell) != 0:
                setattr(l_plus_str,field_l(ell),np.zeros((sum_intersection_dim(ell),1)))
                for c1 in np.arange(1,sum_intersection_dim(ell)+1).reshape(-1):
                    # i = station entrance, j = station exit
                    i = getattr(E_cal_out,(field_E_cal_out(ell)))(c1,1)
                    j = getattr(E_cal_out,(field_E_cal_out(ell)))(c1,2)
                    # Field names of s_s_str, r_s_str and l_str related to the station i,j
                    temp_field_s_s = strcat('s_s_',string(i))
                    temp_field_r_s = strcat('r_s_',string(j))
                    temp_field_l = strcat('l_',string(j))
                    # l_plus calculation
                    getattr[l_plus_str,[temp_field_l]][c3[j]] = getattr(l_str,(temp_field_l))(c3(j)) + T(ell) * (getattr(s_s_str,(temp_field_s_s))(c2(i)) - getattr(r_s_str,(temp_field_r_s))(c3(j)))
                    # Update c2 and c3
                    c2[i] = c2(i) + 1
                    c3[j] = c3(j) + 1
            else:
                setattr(l_plus_str,field_l(ell),[])
        
        ## Cells density
        
        # rho_plus [Nx1] (veh/km):  cells density in the interval k+1
        rho_plus = rho + np.multiply(T / len_,(Phi_plus - Phi_minus))
        
        ## Velocity in the cells
        Phi_avg = (Phi_plus + Phi_minus) / 2
        
        # if the density is equal to 0 then the cars move in free-flow
        ind_rho = find(rho < 0.5)
        v_cell = Phi_avg / rho
        v_cell[ind_rho] = CTM_param.v_bar(ind_rho)
    
    except CTMfunctionError:
        print("Exception in function CTM\n")
        
        ### TODO: PROPER EXCEPTION HANDLING ###

    return rho_plus,l_plus_str,e_plus_str,phi,s,s_s_str,r,r_s_str,Phi_plus,Phi_minus,D,D_on_ramp_s_str,S,v_cell