from OUP import *
import matplotlib.pyplot as plt
import numpy as np
from math import ceil
from sklearn.decomposition import PCA
import pandas as pd

class Visualize:
    
    color_choices = ["blue", "orange", "red", "pink"]

    __use_latex = False
    @staticmethod
    def toggle_latex(val):
        Visualize.__use_latex = val
        if val:
            plt.rcParams['text.usetex'] = True

    @staticmethod
    def sample_point_visualize(samples, legend, title="", file_name="state_plot", center_at_mean=False, fix_axis=False, alpha=1):
        if len(samples) != len(legend):
            raise ValueError("Legend and sample size don't match")
        # TODO: Keep tick values constant
        # TODO: color, label customization

        # Get unique timestep values
        unique_timesteps = sorted(samples[0]['timestep'].unique())

        for sample_points in samples:
            temp = sorted(sample_points["timestep"].unique())
            if temp != unique_timesteps:
                raise ValueError("Different timestamps")

        # Calculate number of subplots based on unique timestep values
        num_subplots = len(unique_timesteps)
        num_cols = min(num_subplots, 3)
        num_rows = int(np.ceil(num_subplots / num_cols))

        # Create figure
        fig = plt.figure(figsize=(20, 6*num_rows))
        axs = []

        x_max, x_min, y_max, y_min, z_max, z_min = float("-inf"), float("inf"), float("-inf"), float("inf"), float("-inf"), float("inf") 
        for j in range(len(samples)):
            sample_points = samples[j]
            # Plot data for each timestep
            for i, timestep in enumerate(unique_timesteps):
                if len(axs) == i:
                    axs.append(fig.add_subplot(num_rows, num_cols, i+1, projection="3d"))
                ax = axs[i]
                data = sample_points[sample_points['timestep'] == timestep]
                if center_at_mean:
                    mean_x = data['x'].mean()
                    mean_y = data['y'].mean()
                    mean_z = data['z'].mean()
                    data.loc[:, "x"] = data["x"].apply(lambda x: x - mean_x)
                    data.loc[:, "y"] = data["y"].apply(lambda y: y - mean_y)
                    data.loc[:, "z"] = data["z"].apply(lambda z: z - mean_z)
                
                ax.scatter(data['x'], data['y'], data['z'], label='timestep={}'.format(timestep), color=Visualize.color_choices[j], alpha=alpha)
                ax.set_xlabel('x (km)')
                ax.set_ylabel('y (km)')
                ax.set_zlabel('z (km)')
                ax.set_title('Timestep {}'.format(timestep))
                ax.legend(legend)
                
                if fix_axis:
                    # Calculate axis limits
                    x_min = min(x_min, data['x'].min())
                    x_max = max(x_max, data['x'].max())

                    y_min = min(y_min, data['y'].min())
                    y_max = max(y_max, data['y'].max())

                    z_min = min(z_min, data['z'].min())
                    z_max = max(z_max, data['z'].max())

                    
                    # Set fixed axis limits


            # Hide unused subplots
            for j in range(i+1, len(axs)):
                axs[j].axis('off')

        if fix_axis:
            for ax in axs:
                ax.set_xlim([x_min, x_max])
                ax.set_ylim([y_min, y_max])
                ax.set_zlim([z_min, z_max])
        # TODO: solve the figure name issue
        plt.suptitle(title, fontweight="bold", fontsize="x-large")
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(use_data_path('visualize/%s.png'%file_name))
        plt.close()

    
    @staticmethod
    def execution_time(time_objs, legend, title="", file_name="execution_time"):
        if len(legend) != len(time_objs):
            raise("Legend size doesnt match given objs")
        
        keys = sorted(time_objs[0]["execution_time_for_step"].keys())
        for i in range(1, len(time_objs)):
            if sorted(time_objs[i]["execution_time_for_step"].keys()) != keys:
                raise("Different timesteps")
        fig, ((ax1), (ax2)) = plt.subplots(2, 1, figsize=(10, 8))

        for i in range(len(time_objs)):
            ax1.plot(keys, [time_objs[i]["execution_time_for_step"][key] for key in keys], color=Visualize.color_choices[i])
        ax1.tick_params(axis="x", labelrotation=90)
        ax1.set_title("Execution time for each timestep")
        ax1.set_xticks(keys)
        ax1.set_xlabel("Time step (s)")
        ax1.set_ylabel("Time taken (s)")

        for i in range(len(time_objs)):
            ax2.barh([i], [time_objs[i]["total_execution_time"]], color=Visualize.color_choices[i]) 
        ax2.invert_yaxis()
        ax2.set_yticks(range(len(legend)), labels=legend)
        ax2.set_xlabel("Time taken (s)")
        ax2.set_title("Total execution time")
        # for bars in ax2.containers:
        #     ax2.bar_label(bars)

        
        fig.legend(legend)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.suptitle(title, fontweight="bold", fontsize="x-large")
        fig.savefig(use_data_path('visualize/%s.png'%(file_name)))  
        plt.close()

    @staticmethod
    def covariance_ellipsoid(dfs, legend, timestep=None, title="", file_name="covariance_ellipsoid", type='wireframe', alpha=1):
        if len(dfs) != len(legend):
            raise "Legend size not same as data size"
        
        timesteps = sorted(dfs[0]["timestep"].unique())
        ntimesteps = len(timesteps)

        for df in dfs:
            temp = sorted(df["timestep"].unique())
            if temp != timesteps:
                raise "Different timesteps"
            
        num_rows = 1 if timestep != None else ceil(ntimesteps/3)
        num_cols = 1 if timestep != None else min(3, ntimesteps)

        fig1 = plt.figure(figsize=(8*num_cols, num_rows*5))
        fig2 = plt.figure(figsize=(8*num_cols, num_rows*5))
        axs1 = []
        axs2 = []

        for j in range(len(dfs)):
            df = dfs[j]
            rellipse = []
            vellipse = []
            mx = None
            my = None
            mz = None
            mvx = None
            mvy = None
            mvz = None
            mix = None
            miy = None
            miz = None
            mivx = None
            mivy = None
            mivz = None
            for i in range(ntimesteps):
                if timestep != None and timesteps[i] != timestep:
                    continue
                df_t = df.loc[df["timestep"] == timesteps[i]]
                cov_r = np.cov(df_t[["x", "y", "z"]].values.T)
                cov_v = np.cov(df_t[["vx", "vy", "vz"]].values.T)
                mean_r = np.mean(df_t[["x", "y", "z"]].values, axis=0)
                mean_v = np.mean(df_t[["vx", "vy", "vz"]].values, axis=0)

                x, y, z = Visualize._get_covarince_ellipsoid(cov_r, mean_r)
                vx, vy, vz = Visualize._get_covarince_ellipsoid(cov_v, mean_v)
                rellipse.append([x, y, z])
                vellipse.append([vx, vy, vz])

                if mx == None:
                    mx = np.max(x)
                    my = np.max(y)
                    mz = np.max(z)
                    mvx = np.max(vx)
                    mvy = np.max(vy)
                    mvz = np.max(vz)
                    mix = np.min(x)
                    miy = np.min(y)
                    miz = np.min(z)
                    mivx = np.min(vx)
                    mivy = np.min(vy)
                    mivz = np.min(vz)
                else:
                    mx = max(mx, np.max(x))
                    my = max(my, np.max(y))
                    mz = max(mz, np.max(z))
                    mvx = max(mvx, np.max(vx))
                    mvy = max(mvy, np.max(vy))
                    mvz = max(mvz, np.max(vz))
                    mix = min(mix, np.min(x))
                    miy = min(miy, np.min(y))
                    miz = min(miz, np.min(z))
                    mivx = min(mivx, np.min(vx))
                    mivy = min(mivy, np.min(vy))
                    mivz = min(mivz, np.min(vz))

            for i in range(len(rellipse)):
                if len(axs1) == i:
                    axs1.append(fig1.add_subplot(num_rows, num_cols, len(axs1)+1, projection="3d"))
                    axs2.append(fig2.add_subplot(num_rows, num_cols, len(axs2)+1, projection="3d"))
                ax1 = axs1[i]
                ax2 = axs2[i]

                ax1.set_title("Position covariance ellipsoid timestep: " + str(timesteps[i]))
                ax2.set_title("Velocity covariance ellipsoid timestep: " + str(timesteps[i]))
                ax1.set_xlabel("x (km)")
                ax1.set_ylabel("y (km)")
                ax1.set_zlabel("z (km)")
                ax2.set_xlabel("x (km/s)")
                ax2.set_ylabel("y (km/s)")
                ax2.set_zlabel("z (km/s)")
                # ax1.set_xlim([mix, mx])
                # ax1.set_ylim([miy, my])
                # ax1.set_zlim([miz, mz])
                # ax2.set_xlim([mivx, mvx])
                # ax2.set_ylim([mivy, mvy])
                # ax2.set_zlim([mivz, mvz])

                if type=='wireframe':    
                    ax1.plot_wireframe(rellipse[i][0], rellipse[i][1], rellipse[i][2], color=Visualize.color_choices[j], alpha=alpha)
                    ax2.plot_wireframe(vellipse[i][0], vellipse[i][1], vellipse[i][2], color=Visualize.color_choices[j], alpha=alpha)
                elif type=='solid':
                    ax1.plot_surface(rellipse[i][0], rellipse[i][1], rellipse[i][2], color=Visualize.color_choices[j], alpha=alpha)
                    ax2.plot_surface(vellipse[i][0], vellipse[i][1], vellipse[i][2], color=Visualize.color_choices[j], alpha=alpha)
                else:
                    raise ValueError("type should be either 'solid' or 'wireframe'")
                
                ax1.legend(legend)
                ax2.legend(legend)


        fig1.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig1.suptitle(title, fontweight="bold", fontsize="x-large")
        fig1.savefig(use_data_path("visualize/%s_position.png"%file_name))
        fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig2.suptitle(title, fontweight="bold", fontsize="x-large")
        fig2.savefig(use_data_path("visualize/%s_velocity.png"%file_name))

        plt.close(fig1)
        plt.close(fig2)

    @staticmethod
    def _get_covarince_ellipsoid(cov, mean, nstd=3):
        assert cov.shape==(3,3)

        # Find and sort eigenvalues to correspond to the covariance matrix
        eigvals, eigvecs = np.linalg.eigh(cov)
        idx = np.sum(cov,axis=0).argsort()
        eigvals_temp = eigvals[idx]
        idx = eigvals_temp.argsort()
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:,idx]

        # Set of all spherical angles to draw our ellipsoid
        n_points = 20
        theta = np.linspace(0, 2*np.pi, n_points)
        phi = np.linspace(0, np.pi, n_points)

        # Width, height and depth of ellipsoid
        rx, ry, rz = nstd * np.sqrt(eigvals)

        # Get the xyz points for plotting
        # Cartesian coordinates that correspond to the spherical angles:
        X = rx * np.outer(np.cos(theta), np.sin(phi))
        Y = ry * np.outer(np.sin(theta), np.sin(phi))
        Z = rz * np.outer(np.ones_like(theta), np.cos(phi))

        # Rotate ellipsoid for off axis alignment
        old_shape = X.shape
        # Flatten to vectorise rotation
        X,Y,Z = X.flatten(), Y.flatten(), Z.flatten()
        X,Y,Z = np.matmul(eigvecs, np.array([X,Y,Z]))
        X,Y,Z = X.reshape(old_shape), Y.reshape(old_shape), Z.reshape(old_shape)
    
        # Add in offsets for the mean
        # X = X + mean[0]
        # Y = Y + mean[1]
        # Z = Z + mean[2]
        
        return X, Y, Z
    
    @staticmethod
    def covariance_and_samples(df, timestep=None, title="", file_name=""):
        timesteps = sorted(df["timestep"].unique())
        ntimesteps = len(timesteps)

        fig = plt.figure(figsize=(15, ntimesteps*5 if timestep==None else 5))
        axs = []

        rellipse = []
        vellipse = []
        mx = None
        my = None
        mz = None
        mvx = None
        mvy = None
        mvz = None
        mix = None
        miy = None
        miz = None
        mivx = None
        mivy = None
        mivz = None
        for i in range(ntimesteps):
            if timestep != None and timesteps[i] != timestep:
                continue
            df_t = df.loc[df["timestep"] == timesteps[i]]
            cov_r = np.cov(df_t[["x", "y", "z"]].values.T)
            cov_v = np.cov(df_t[["vx", "vy", "vz"]].values.T)
            mean_r = np.mean(df_t[["x", "y", "z"]].values, axis=0)
            mean_v = np.mean(df_t[["vx", "vy", "vz"]].values, axis=0)

            x, y, z = Visualize._get_covarince_ellipsoid(cov_r, mean_r)
            vx, vy, vz = Visualize._get_covarince_ellipsoid(cov_v, mean_v)
            rellipse.append([x, y, z, timesteps[i]])
            vellipse.append([vx, vy, vz, timesteps[i]])

            if mx == None:
                mx = np.max(x)
                my = np.max(y)
                mz = np.max(z)
                mvx = np.max(vx)
                mvy = np.max(vy)
                mvz = np.max(vz)
                mix = np.min(x)
                miy = np.min(y)
                miz = np.min(z)
                mivx = np.min(vx)
                mivy = np.min(vy)
                mivz = np.min(vz)
            else:
                mx = max(mx, np.max(x))
                my = max(my, np.max(y))
                mz = max(mz, np.max(z))
                mvx = max(mvx, np.max(vx))
                mvy = max(mvy, np.max(vy))
                mvz = max(mvz, np.max(vz))
                mix = min(mix, np.min(x))
                miy = min(miy, np.min(y))
                miz = min(miz, np.min(z))
                mivx = min(mivx, np.min(vx))
                mivy = min(mivy, np.min(vy))
                mivz = min(mivz, np.min(vz))

        for i in range(len(rellipse)):
            if len(axs) == i:
                axs.append([fig.add_subplot(len(rellipse), 2, 2*i+1, projection="3d"), fig.add_subplot(len(rellipse), 2, 2*i+2, projection="3d")])
            ax1 = axs[i][0]
            ax2 = axs[i][1]

            if timestep == None:
                ax1.set_title("Position timestep: " + str(timesteps[i]))
                ax2.set_title("Velocity timestep: " + str(timesteps[i]))
            else:
                ax1.set_title("Position")
                ax2.set_title("Velocity")

            ax1.set_xlabel("x (km)")
            ax1.set_ylabel("y (km)")
            ax1.set_zlabel("z (km)")
            ax2.set_xlabel("x (km/s)")
            ax2.set_ylabel("y (km/s)")
            ax2.set_zlabel("z (km/s)")

            # ax1.set_xlim([mix, mx])
            # ax1.set_ylim([miy, my])
            # ax1.set_zlim([miz, mz])
            # ax2.set_xlim([mivx, mvx])
            # ax2.set_ylim([mivy, mvy])
            # ax2.set_zlim([mivz, mvz])

            data = df[df["timestep"] == rellipse[i][3]]
            ax1.plot_wireframe(rellipse[i][0], rellipse[i][1], rellipse[i][2], color=Visualize.color_choices[0], zorder=1)
            ax1.scatter(data['x'], data['y'], data['z'], color=Visualize.color_choices[1], zorder=2)
            ax1.set_xlabel('x')
            ax1.set_ylabel('y')
            ax1.set_zlabel('z')
            ax1.set_title('Timestep {}'.format(rellipse[i][3]))
            
            ax2.plot_wireframe(vellipse[i][0], vellipse[i][1], vellipse[i][2], color=Visualize.color_choices[0], zorder=1)
            ax2.scatter(data['vx'], data['vy'], data['vz'], color=Visualize.color_choices[1], zorder=2)
            ax2.set_xlabel('x')
            ax2.set_ylabel('y')
            ax2.set_zlabel('z')
            ax2.set_title('Timestep {}'.format(rellipse[i][3]))
        


        fig.tight_layout()
        fig.suptitle(title)
        fig.savefig(use_data_path("visualize/%s.png"%file_name))
        plt.close()

    @staticmethod
    def comparision_table(samples, legend, title1="", title2="", file_name="comparision_table"):
        if len(samples) != len(legend):
            raise ValueError("Legend and sample size don't match")
        # TODO: Keep tick values constant
        # TODO: color, label customization

        # Get unique timestep values
        unique_timesteps = sorted(samples[0]['timestep'].unique())

        for sample_points in samples:
            temp = sorted(sample_points["timestep"].unique())
            if temp != unique_timesteps:
                raise ValueError("Different timestamps")

        '''
        For each timestep, generate a table, structure it as follows
        S.No|Parameter           | MC    | ESPT  | AESPT 
        1.  |mean sample pt      |       |       |
        2.  |distance b/w mean   |       |       |
            |point MC            |       |       |
        3.  |distance b/w mean   |       |       |
            |point ESPT          |       |       |
        4.  |distance b/w mean   |       |       |
            |point AESPT         |       |       |
        5.  |standard deviation  |       |       |
            |of sample points    |       |       |
        '''

        table = {}

        for i, timestep in enumerate(unique_timesteps):
            table_for_timestep = {}
            table_for_timestep['r_mean_sample_pts'] = []
            table_for_timestep['r_std_dev_sample_pts'] = []
            table_for_timestep['v_mean_sample_pts'] = []
            table_for_timestep['v_std_dev_sample_pts'] = []
            for j in range(len(samples)):
                sample_points = samples[j]
                data = sample_points[sample_points['timestep'] == timestep]
                table_for_timestep['r_mean_sample_pts'].append(data[["x", "y", "z"]].mean())
                table_for_timestep['r_std_dev_sample_pts'].append(data[["x", "y", "z"]].std())
                table_for_timestep['v_mean_sample_pts'].append(data[["vx", "vy", "vz"]].mean())
                table_for_timestep['v_std_dev_sample_pts'].append(data[["vx", "vy", "vz"]].std())
        
            table[timestep] = table_for_timestep

        # for i, timestep in enumerate(unique_timesteps):
        #     for j in range(len(samples)):
        #         print()
        #         print(legend[j], ' -> ', timestep)
        #         print('Mean sample point: ', table[timestep]['mean_sample_pts'][j])
        #         print('Standard deviation of the sample point: ', table[timestep]['std_dev_sample_pts'][j])
        
        fig1, ax1 = plt.subplots(2, len(samples)-1, figsize=(8*(len(samples)-1), 10))
        for i in range(len(samples)-1):
            axtr = ax1[0][i]
            axtv = ax1[1][i]

            bars1 = axtr.bar(range(len(unique_timesteps)), [np.linalg.norm(table[t]["r_mean_sample_pts"][0] - table[t]["r_mean_sample_pts"][i+1]) for t in unique_timesteps])
            # ax1.bar_label(bars1)
            axtr.set_ylabel("Distance between mean position(km)")
            axtr.set_xlabel("Timestep (s)")
            axtr.set_xticks(range(len(unique_timesteps)), labels=unique_timesteps)
            axtr.set_title("Mean position distance over time: %svs%s"%(legend[0], legend[i+1]))

            bars2 = axtv.bar(range(len(unique_timesteps)), [np.linalg.norm(table[t]["v_mean_sample_pts"][0] - table[t]["v_mean_sample_pts"][i+1]) for t in unique_timesteps])
            # ax1.bar_label(bars1)
            axtv.set_ylabel("Difference between mean velocity(km/s)")
            axtv.set_xlabel("Timestep (s)")
            axtv.set_xticks(range(len(unique_timesteps)), labels=unique_timesteps)
            axtv.set_title("Mean velocity diffference over time: %svs%s"%(legend[0], legend[i+1]))
        

            
        fig2, ax2 = plt.subplots(2, 1, figsize=(10, 10))
        for i in range(len(samples)):
            ax2[0].plot(unique_timesteps, [np.linalg.norm(table[ts]["r_std_dev_sample_pts"][i]) for ts in unique_timesteps], color=Visualize.color_choices[i])
        ax2[0].set_xticks(unique_timesteps)
        ax2[0].set_xlabel("Timestep (s)")
        if Visualize.__use_latex:
            ax2[0].set_ylabel(r'\[\sqrt{\sigma_x^2 + \sigma_y^2 + \sigma_z^2}\] (km)')
        else:
            ax2[0].set_ylabel('Norm of std dev of position (km)')
        ax2[0].legend(legend)
        ax2[0].set_title("Standard deviation evolution of position with propagation")

        for i in range(len(samples)):
            ax2[1].plot(unique_timesteps, [np.linalg.norm(table[ts]["v_std_dev_sample_pts"][i]) for ts in unique_timesteps], color=Visualize.color_choices[i])
        ax2[1].set_xticks(unique_timesteps)
        ax2[1].set_xlabel("Timestep (s)")
        if Visualize.__use_latex:
            ax2[1].set_ylabel(r'\[\sqrt{\sigma_x^2 + \sigma_y^2 + \sigma_z^2}\] (km/s)')
        else:
            ax2[1].set_ylabel('Norm of std dev of velocity (km/s)')
        ax2[1].legend(legend)
        ax2[1].set_title("Standard deviation evolution of velocity with propagation")

        fig1.suptitle(title1, fontweight="bold", fontsize="x-large")
        fig2.suptitle(title2, fontweight="bold", fontsize="x-large")

        fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig2.savefig(use_data_path("visualize/%s_stddev.png"%file_name))

        fig1.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig1.savefig(use_data_path("visualize/%s_mean.png"%file_name))

        plt.close(fig1)
        plt.close(fig2)
        
            
        for i, timestep in enumerate(unique_timesteps):
                print()
                print('MC', ' -> ', timestep)
                print('Norm of mean_MC - mean_ESPT: ',np.linalg.norm(table[timestep]['r_mean_sample_pts'][0] - table[timestep]['r_mean_sample_pts'][1]))
                print('Norm of mean_MC - mean_AESPT: ',np.linalg.norm(table[timestep]['r_mean_sample_pts'][0] - table[timestep]['r_mean_sample_pts'][2]))
                print()
                print('ESPT', ' -> ', timestep)
                print('Norm of mean_ESPT - mean_MC: ',np.linalg.norm(table[timestep]['r_mean_sample_pts'][1] - table[timestep]['r_mean_sample_pts'][0]))
                print('Norm of mean_ESPT - mean_AESPT: ',np.linalg.norm(table[timestep]['r_mean_sample_pts'][1] - table[timestep]['r_mean_sample_pts'][2]))
                print()
                print('AESPT', ' -> ', timestep)
                print('Norm of mean_AESPT - mean_MC: ',np.linalg.norm(table[timestep]['r_mean_sample_pts'][2] - table[timestep]['r_mean_sample_pts'][0]))
                print('Norm of mean_AESPT - mean_ESPT: ',np.linalg.norm(table[timestep]['r_mean_sample_pts'][2] - table[timestep]['r_mean_sample_pts'][1]))

    
    @staticmethod
    def PCA(df, title="pca"):
        # You must normalize the data before applying the fit method
        df_normalized=(df - df.mean()) / df.std()
        pca = PCA(n_components=df.shape[1])
        pca.fit(df_normalized)

        # Reformat and view results
        loadings = pd.DataFrame(pca.components_.T,
        columns=['PC%s' % _ for _ in range(len(df_normalized.columns))],
        index=df.columns)
        print(loadings)

        fig, ax = plt.subplots()
        print(pca.explained_variance_ratio_)
        ax.bar([0, 1, 2], pca.explained_variance_ratio_)
        ax.set_xticks([0, 1, 2], labels=["x", "y", "z"])
        ax.set_ylabel('Explained Variance')
        ax.set_xlabel('Components')
        fig.savefig(use_data_path("visualize/%s.png"%title))
        fig.suptitle("Initial PCA")