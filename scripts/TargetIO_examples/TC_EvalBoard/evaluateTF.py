import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pdb

#outdir = "/home/cta/luigi/TC_EvalBoard/data/20160307/1VpedDAC_SBias1050/"
#infile = outdir+"TFdata.npz"

def DAC2mV(DAC_array):
    """
    based on measure performed by LT using TC eval board
    :param: array of DAC Vped values
    :return: array of Vped values in mV
    """
    mV_array = -2.76 + 0.611 * DAC_array
    return mV_array


class TargetTF:
    def __init__(self, infile):
        data = np.load(infile)
        try:
            self.tf_array = data['arr_0']
            if np.any(self.tf_array == 0):
                print "WARNING: the TF array contains 0 values"
            Vped_DAC = data['arr_1']
            self.Vped = DAC2mV(Vped_DAC)
            if not len(self.Vped) == len(self.tf_array):
                print "mismatch between array size in input file"
            self.nvped, self.nsample, self.nevent = np.shape(self.tf_array)
            self.tfs = np.average(self.tf_array, axis=2)
        except:
            print "failed to load TF data from input file"

    def plot_TF(self, outdir,xmin=0,xmax=-1):
        fig = plt.figure("Transfer Function")
        ax = plt.subplot(111)
        ax.set_xlabel("Vped (mV)", fontsize=14)
        ax.set_ylabel("ADC counts", fontsize=14)
        for cell in range(self.nsample):
            ax.plot(self.Vped, self.tfs[:, cell])
        plt.axvline(self.Vped[xmin],linestyle='--',color='b')
        plt.axvline(np.average(self.Vped[xmax]),linestyle='--',color='b')
        fig.savefig(outdir + "TF.png")
        fig.clf()

    def plot_DC_noise(self, outdir,noise_cut=4, xmin=0,xmax=-1,npoint=2):
        fig = plt.figure("DC noise")
        ax = plt.subplot(111)
        ax.set_xlabel("Vped (mV)", fontsize=14)
        #ax.set_ylabel("RMS (mV)", fontsize=14)
        ax.set_ylabel("RMS (ADC counts)", fontsize=14)
        #ax.set_xlim(600, 2500)
        ax.set_ylim(0, 6)
        rmsADC = np.std(self.tf_array, axis=2)
        slopes = [[linregress(self.Vped[vped - npoint:vped + npoint+1], self.tfs[vped - npoint:vped + npoint+1, cell])[0] for cell in
                   range(self.nsample)] for vped in range(npoint, self.nvped - npoint)]
        slopes = np.array(slopes)
        slopes[slopes<1.e-6] = 1.e-6
        print np.average(slopes[slopes>1.e-6])
        rmsmV = rmsADC[npoint:-npoint, :] / slopes
        for cell in range(self.nsample):
            ax.plot(self.Vped, rmsADC)
            #ax.plot(self.Vped[npoint:-npoint], rmsmV[:, cell])
        plt.axvline(self.Vped[xmin],linestyle='--',color='k')
        plt.axvline(np.average(self.Vped[xmax]),linestyle='--',color='k')
        #evaluate new range
        avg_noise = np.average(rmsmV,axis=1)
        #pdb.set_trace()
        #print self.Vped[npoint:-npoint]
        print avg_noise
        goodrange = np.where(avg_noise<noise_cut)
        x1 = goodrange[0][0]+npoint
        print "x1", x1
        xmin = np.maximum(xmin,x1)
        #print "xmin", xmin
        x2 = goodrange[0][-1]+npoint
        #print "x2", x2
        xmax = np.minimum(xmax,x2)
        #print "xmax", xmax
        plt.axvline(self.Vped[xmin-1],linestyle='--',color='b')
        plt.axvline(np.average(self.Vped[xmax+1]),linestyle='--',color='b')

        fig.savefig(outdir + "DCnoise.png")
        fig.clf()
        return xmin, xmax,rmsmV

    def plot_nonlinearity(self,outdir, bound=10, lindev_cut=400):

        #first pass, to evaluate dynamic range span
        lininterp = np.zeros(np.shape(self.tfs))
        for cell in range(self.nsample):
            slope, intercept, r_value, p_value, std_err  = linregress(self.Vped[bound:-bound], self.tfs[bound:-bound, cell])
            lininterp[:, cell] = intercept + slope * self.Vped
        lindev = self.tfs - lininterp
        xmin = np.min([np.min(np.where(np.abs(lindev[:,cell])<lindev_cut)) for cell in range(self.nsample)])
        xmax = [np.max(np.where(np.abs(lindev[:,cell])<lindev_cut)) for cell in range(self.nsample)]
        x2 = [np.where(self.tfs[:,cell]==np.max(self.tfs[:,cell]))[0][0] for cell in range(self.nsample)]
        xmax = [np.maximum(xmax[cell],x2[cell]) for cell in range(self.nsample)]
        xmax = np.array(xmax)

        #second pass, now do linear fit and calculate non linearity over range above
        lininterp = np.zeros(np.shape(self.tfs))
        for cell in range(self.nsample):
            slope, intercept, r_value, p_value, std_err  = linregress(self.Vped[xmin:np.min(xmax)], self.tfs[xmin:np.min(xmax), cell])
            lininterp[:, cell] = intercept + slope * self.Vped
        lindev = self.tfs - lininterp

        fig = plt.figure("Nonlinearity")
        ax = plt.subplot(111)
        ax.set_xlabel("Vped (mV)", fontsize=14)
        ax.set_ylabel("Deviation from linear fit (ADC counts)", fontsize=14)
        ax.set_ylim(-500,500)
        for cell in range(self.nsample):
            ax.plot(self.Vped, lindev[:, cell])
        plt.axvline(self.Vped[xmin-1],linestyle='--',color='k')
        plt.axvline(np.average(self.Vped[np.minimum(xmax+1,len(self.Vped)-1)]),linestyle='--',color='k')
        fig.savefig(outdir + "Nonlinearity.png")
        fig.clf()

        INL = np.array([np.max(lindev[xmin:np.min(xmax),cell]) for cell in range(self.nsample)])
        INL = np.average(INL)

        return xmin, xmax, INL
    
    def evaluate(self,outdir,npoint=2):
        xmin, xmax, INL = self.plot_nonlinearity(outdir)
        xmin, xmax, DC_noise = self.plot_DC_noise(outdir,xmin=xmin,xmax=xmax,npoint=npoint)
        self.plot_TF(outdir,xmin=xmin,xmax=xmax)
        print "Vmin, Vmax: ", self.Vped[xmin-1],  np.average([self.Vped[xmax[cell]+1] for cell in range(self.nsample)])
        drange = np.average([self.Vped[xmax[cell]+1]-self.Vped[xmin-1] for cell in range(self.nsample)])
        print "dynamic range", drange, "mV"
        print "INL", INL, "ADC counts"
        avdcnoise = np.average(DC_noise[xmin-npoint:np.min(xmax)-npoint,:])
        print "DC noise", avdcnoise, "mV"
        print "effective dynamic range", np.log2(drange/avdcnoise), "bits"
        return drange, INL, avdcnoise, np.log2(drange/avdcnoise)



#tf = TargetTF(infile)
#drange, INL, avdcnoise, eff_drange = tf.evaluate(outdir,npoint=10)

