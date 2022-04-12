    # print(sigma_ksi[:,0])
    # quit()

    # S_ksi, T_ksi = sp.linalg.eig(P_ksi)
    # S_ksi = np.real(S_ksi)

    # for i in range(n_ksi):
    #     mat = np.zeros((n_ksi,n_ksi))
    #     mat[i,i] = np.sqrt(S_ksi[i])
    #     sigma_ksi[:,2*i+0] = ksi + T_ksi.T @ mat
    #     sigma_ksi[:,2*i+1] = ksi - T_ksi.T @ mat

    # ## P_m term 
    # Tp, Sp = sp.linalg.eig(P_m)
    # ## sig_ww term
    # Tw, Sw = sp.linalg.eig(sig_ww)

    # Sp = np.real(Sp)

    # sigma_ksi[:,1] = ksi + Tp.T @ np.diag(np.sqrt(Sp[0]), 0, 0) * n_ksi
    # sigma_ksi[:,2] = ksi - Tp.T @ np.diag(np.sqrt(Sp[0]), 0, 0) * n_ksi
    # sigma_ksi[:,3] = ksi + Tp.T @ np.diag(0, np.sqrt(Sp[1]), 0) * n_ksi
    # sigma_ksi[:,4] = ksi - Tp.T @ np.diag(0, np.sqrt(Sp[1]), 0) * n_ksi 
    # sigma_ksi[:,5] = ksi + Tp.T @ np.diag(0, 0, np.sqrt(Sp[2])) * n_ksi
    # sigma_ksi[:,6] = ksi - Tp.T @ np.diag(0, 0, np.sqrt(Sp[2])) * n_ksi

    # for i in range (4):
    #     sigma_ksi[:,7+2*i+0] = ksi + np.sqrt(sig_vv[i]) * n_ksi
    #     sigma_ksi[:,7+2*i+1] = ksi - np.sqrt(sig_vv[i]) * n_ksi
    
    # sigma_ksi[:,15] = ksi + Tw.T @ np.diag(np.sqrt(Sw[0]), 0) * n_ksi
    # sigma_ksi[:,16] = ksi - Tw.T @ np.diag(np.sqrt(Sw[0]), 0) * n_ksi
    # sigma_ksi[:,17] = ksi + Tw.T @ np.diag(0, np.sqrt(Sw[1])) * n_ksi
    # sigma_ksi[:,18] = ksi - Tw.T @ np.diag(0, np.sqrt(Sw[1])) * n_ksi