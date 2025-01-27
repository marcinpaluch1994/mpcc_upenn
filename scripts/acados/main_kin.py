import numpy as np
import yaml
from matplotlib import pyplot as plt

import Bezier
from acados_settings import acados_settings_kin
from python_sim_utils import plotter, compute_objective


def main_kin():
    np.set_printoptions(precision=3)

    # model parameters
    paramfile = "modelparams.yaml"
    #load global constant model parameters
    with open(paramfile) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    lencar = lf + lr

    #sim parameters
    Tsim = 4
    Tf = 1
    N = 20
    Qc = 0.1
    Ql = 100
    Q_theta = 1
    R_d = 0.1
    R_delta = 0.1
    r = 0.2 #trackwidth
    Nsim = int(np.floor(N/Tf*Tsim))


    track_lu_table, smax = Bezier.generatelookuptable("tracks/simpleoval")
    trk_plt = plotter(track_lu_table, smax, r, lencar)
    trk_plt.plot_track()
    plt.show()

    constraints, model, acados_solver, ocp = acados_settings_kin(Tf, N,  paramfile)

    #starting position in track startidx = theta0[m] * 100 [pts/m]
    startidx = 1050

    trackvars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    xvars = ['posx', 'posy', 'phi', 'vx', 'theta', 'd', 'delta']
    uvars = ['ddot', 'deltadot', 'thetadot']
    car_soln = []
    xt0 = track_lu_table[startidx,trackvars.index('xtrack')]
    yt0 = track_lu_table[startidx,trackvars.index('ytrack')]
    phit0 = track_lu_table[startidx,trackvars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx,trackvars.index('sval')]
    x0 = np.array([xt0, yt0, phit0, 1, theta_hat0, 0, 0])

    ############################################################################
    #initialization for theta values
    iter = 30
    theta_old = theta_hat0*np.ones((N,))
    x_current = np.tile(x0,(N,1))

    for idx in range(iter):
        #get track linearization
        index_lin_points = 100 * theta_old
        index_lin_points = index_lin_points.astype(np.int32)
        track_lin_points = track_lu_table[index_lin_points,:]

        for stageidx in range(N):
            p_val = np.array([track_lin_points[stageidx,trackvars.index('xtrack')],
                                track_lin_points[stageidx,trackvars.index('ytrack')],
                                track_lin_points[stageidx,trackvars.index('phitrack')],
                                track_lin_points[stageidx,trackvars.index('sin(phi)')],
                                track_lin_points[stageidx,trackvars.index('cos(phi)')],
                                track_lin_points[stageidx,trackvars.index('sval')],  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta,
                                r-0.5*lencar
                                ])
            #print("stage idx: ",stageidx,"pval: ",p_val[:-6])
            x_val = x_current[stageidx]
            acados_solver.set(stageidx,"p", p_val)
            acados_solver.set(stageidx,"x", x_val)
        #constrain x0
        acados_solver.set(0, "lbx", x0)
        acados_solver.set(0, "ubx", x0)

        #solve problem
        status = acados_solver.solve()
        #extract theta values
        for idx_sol in range(N):
            xsol = acados_solver.get(idx_sol,"x")
            x_current[idx_sol,:] = xsol

        theta_current = x_current[:,4]

        #compute difference
        theta_diff = np.sum(np.abs(theta_current-theta_old))
        print("theta init difference: ", theta_diff)
        print("theta values", theta_current)
        theta_old = theta_current
    ############################################################################
    #setup using estimated initial trajectory
    theta_vals = theta_current

    step_sol_x_arr = x_current

    print("plotting initialization")
    trk_plt.plot_horizon(theta_vals, step_sol_x_arr[:, :3])
    plt.pause(0.1)
    input("hit [enter] to continue.")
    plt.pause(0.1)
    trk_plt.clear_horizion()


    #list storing visited states
    x0vals = []
    laps = 0
    ##########################SIMULATION#######################################
    for simidx in range(Nsim):

        step_sol_x = []
        step_sol_u = []

        theta_old = theta_vals
        #get track linearization
        index_lin_points = 100 * theta_old - 100*laps*smax
        index_lin_points = index_lin_points.astype(np.int32)
        print("track linearized around entries:", index_lin_points )
        track_lin_points = track_lu_table[index_lin_points,:]

        #######################################################################
        #set params and warmstart
        for stageidx in range(N-1):
            p_val = np.array([track_lin_points[stageidx,trackvars.index('xtrack')],
                                track_lin_points[stageidx,trackvars.index('ytrack')],
                                track_lin_points[stageidx,trackvars.index('phitrack')],
                                track_lin_points[stageidx,trackvars.index('sin(phi)')],
                                track_lin_points[stageidx,trackvars.index('cos(phi)')],
                                track_lin_points[stageidx,trackvars.index('sval')] + laps*smax,  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta,
                                r-0.5*lencar
                                ])
            #print("stage idx: ",stageidx)
            #print("pval: ",p_val[:-5])

            acados_solver.set(stageidx,"p", p_val)
            acados_solver.set(stageidx, "x", step_sol_x_arr[stageidx+1])

        #last stage copy old solution for init
        stageidx = N-1
        p_val = np.array([track_lin_points[stageidx,trackvars.index('xtrack')],
                            track_lin_points[stageidx,trackvars.index('ytrack')],
                            track_lin_points[stageidx,trackvars.index('phitrack')],
                            track_lin_points[stageidx,trackvars.index('sin(phi)')],
                            track_lin_points[stageidx,trackvars.index('cos(phi)')],
                            track_lin_points[stageidx,trackvars.index('sval')] + laps*smax,  #aka theta_hat
                            Qc,
                            Ql,
                            Q_theta,
                            R_d,
                            R_delta,
                            r-0.5*lencar
                            ])
        #print("stage idx: ",stageidx,"pval: ",p_val[:-6])
        acados_solver.set(stageidx,"p", p_val)
        acados_solver.set(stageidx, "x", step_sol_x_arr[stageidx])
        #######################################################################
        #constrain x0
        acados_solver.set(0, "lbx", x0)
        acados_solver.set(0, "ubx", x0)

        status = acados_solver.solve()
        #acados_solver.print_statistics()
        #print(status)

        #extract solution
        x0 = acados_solver.get(1,"x")
        x0vals.append(x0)

        for idx_sol in range(N):
            xsol = acados_solver.get(idx_sol,"x")
            #print("stage:",idx_sol," xsol:", xsol)
            step_sol_x.append(xsol)
            step_sol_u.append(acados_solver.get(idx_sol,"u"))

        step_sol_x_arr = np.array(step_sol_x)
        step_sol_u_arr = np.array(step_sol_u)

        theta_vals = step_sol_x_arr[:, 4]
        print("theta vals", theta_vals)

        objective = compute_objective(Tf/float(N),
                    Qc,
                    Ql,
                    Q_theta,
                    R_d,
                    R_delta,
                    theta_vals,
                    theta_old,
                    step_sol_x_arr[:, :2],
                    step_sol_u_arr,
                    track_lin_points[:,trackvars.index('xtrack'):trackvars.index('ytrack')+1],
                    track_lin_points[:,trackvars.index('phitrack')]
                    )


        #plotting result

        print("objective value", objective)
        trk_plt.plot_horizon(theta_vals, step_sol_x_arr[:, :3])
        trk_plt.plot_input_state_traj(step_sol_x_arr, step_sol_u_arr, xvars, uvars)
        #plt.show()

        plt.pause(0.1)
        # input("hit [enter] to continue.")
        plt.pause(0.1)
        trk_plt.clear_horizion()
        trk_plt.clear_input_state_traj()

        #preparation for next timestep
        theta_vals = np.hstack((step_sol_x_arr[1:, 4], step_sol_x_arr[-1, 4]+0.1))

        if theta_vals[0] > (laps+1)*smax :
            print("#################################RESET###############################")
            laps = laps + 1


    ###############################/SIMULATION##################################
    trk_plt.plot_traj(np.array(x0vals))
    plt.show()
    #np.savetxt("full_sol_x_log.csv", )
        #print(simidx)
    return 0
