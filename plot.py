import numpy as np
import matplotlib.pyplot as plt

def plot_serial():
    nside=128
    T=np.loadtxt("T_out.txt")

    plt.imshow(T)
    plt.title("T for nside="+str(nside))
    plt.colorbar()
    plt.show()
    plt.close()

def plot_mpi():
    nside=[128,256]#,512]
    nproc=[1,2,4,8,16]
    xymax=np.pi

    for ns in nside:
        T_grid=np.zeros((ns,ns))
        for npr in nproc:
            for i in range(npr):
                data=np.loadtxt("OutputDatafiles/T_out_nside"+str(ns)+"_process"+str(i)+"_numproc"+str(npr)+".txt")
                nx=data.shape[0]-2
                i_min=nx*i
                i_max=nx*(i+1)
                data=data[1:data.shape[0]-1,:]
                T_grid[i_min:i_max, :]=data

            T_grid=np.transpose(T_grid) #so that x plots on horizontal axis

            plt.imshow(T_grid, origin='lower', extent=[0,xymax,0,xymax])#, cmap='RdBu')
            plt.colorbar()
            plt.title("T for "+str(npr)+" processor(s), nside="+str(ns))
            plt.savefig("Plots/T_"+str(npr)+"processors_nside"+str(ns)+".png")
            plt.xlabel('x')
            plt.ylabel('y')
            plt.show()
            plt.close()





if __name__=="__main__":
    #plot_serial()
    #plot_omp()
    plot_mpi()
