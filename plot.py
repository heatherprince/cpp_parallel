import numpy as np
import matplotlib.pyplot as plt

def plot_T(T_grid, title, saveas):
    xymax=np.pi
    plt.imshow(T_grid, origin='lower', extent=[0,xymax,0,xymax])#, cmap='RdBu')
    plt.colorbar()
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(saveas)
    plt.show()
    plt.close()


def plot_serial(nside):
    for ns in nside:
        T=np.loadtxt("OutputDatafiles/T_out_nside"+str(ns)+"_serial.txt")
        T_grid=np.transpose(T)
        plot_T(T_grid, "Serial: T calculated for nside="+str(ns), "Plots/T_serial_nside"+str(ns)+".png")

def plot_mpi(nside, nproc):
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
            plot_T(T_grid, "MPI: T for "+str(npr)+" processor(s), nside="+str(ns), "Plots/T_"+str(npr)+"processors_nside"+str(ns)+".png")


def plot_omp(nside, nthreads):
    for ns in nside:
        T_grid=np.zeros((ns,ns))
        for nth in nthreads:
            T=np.loadtxt("OutputDatafiles/T_out_nside"+str(ns)+"_omp_nthreads"+str(nth)+".txt")
            T_grid=np.transpose(T) #so that x plots on horizontal axis
            plot_T(T_grid, "OpenMP: T for "+str(nth)+" thread(s), nside="+str(ns), "Plots/T_"+str(nth)+"threads_nside"+str(ns)+".png")




if __name__=="__main__":
    plot_serial([128])
    plot_omp([128], [1,4])
    plot_mpi([128], [1,4])
