import numpy as np

tol = 10 ** - 5
vel_diff = (1.0 / vel_cell_st(:,k) - 1.0 / CTM_param.v_bar)
ind_vel_diff = find(vel_diff < tol)
vel_diff[ind_vel_diff] = np.zeros((ind_vel_diff.shape,ind_vel_diff.shape))
Delta_st[:,k] = np.multiply(CTM_param.len,vel_diff)