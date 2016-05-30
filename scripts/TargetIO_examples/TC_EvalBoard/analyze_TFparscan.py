import target_io
import target_driver
import numpy as np
import matplotlib.pyplot as plt
import evaluateTF
import gc

datadir1="/home/cta/luigi/TC_EvalBoard/data/20160405/test/"
#parname = 'SBias'
#parvals = np.arange(900,1200,10)
#parname = "RampSignalDuration"
#parvals = np.arange(350,800,10)
#parname = "Isel"
#parvals = np.arange(2200,2500,400)
#parvals = np.arange(1800,2230,10)
#parvals = np.arange(2010,2500,500)
# parname = "PUbias"
# parvals = np.arange(3008,3018,10)
# parname = "Isel"
# parvals = np.arange(2170,2180,10)
#parname = "PUbias"
#parvals = np.arange(3050,3060,10)
#parvals = np.arange(2950,3410,10)
#parname = "CMPbias2"
#parvals = np.arange(425,835,10)
parname = "Isel"
parvals = np.arange(2400,2404,4)

Vped_DAC = np.arange(495,4096,120)
nsample=96
nevent=300

drange = np.zeros(len(parvals))
INL = np.zeros(len(parvals))
avdcnoise = np.zeros(len(parvals))
eff_drange = np.zeros(len(parvals))

outfile = open(datadir1+"results.txt","w")
outfile.write(parname+" dynamic range (mV) INL (ADC counts) DC noise (mV) efective dynamic range (bits) \n")

for ipar, parval in enumerate(parvals):
    datadir = datadir1+'{}_{:04d}/'.format(parname,parval)
    print "processing data in", datadir
    tf_array=np.zeros([len(Vped_DAC),nsample,nevent])
    #extract values and save to np array
    for s, Vped in enumerate(Vped_DAC):
        ampl = np.zeros([nsample,nevent])
        #print "processing TF data for Vped {} DAC counts".format(Vped)
        filename=datadir+"TFdata_VpedDAC_{:04d}.fits".format(Vped)
        reader = target_io.EventFileReader(filename)
        #print "done"
        for ievt in range(nevent):
            rawdata = reader.GetEventPacket(ievt,0)#second entry is asic/channel, only one enabled in this case
            packet = target_driver.DataPacket()
            packet.Assign(rawdata, reader.GetPacketSize())
            wf = packet.GetWaveform(0)
            for sample in range(nsample):
                ampl[sample,ievt] = wf.GetADC(sample+12)
        tf_array[s,:,:]= ampl[:,:]
        reader.Close()
    np.savez_compressed(datadir+"TFdata.npz",tf_array,Vped_DAC)
    #evaluate TF
    tf = evaluateTF.TargetTF(datadir+"TFdata.npz")
    tf.plot_TF(outdir=datadir)
    #drange1, INL1, avdcnoise1, eff_drange1 = tf.evaluate(datadir,npoint=2)
    #drange[ipar]  = drange1
    #INL[ipar] = INL1
    #avdcnoise[ipar] = avdcnoise1
    #eff_drange[ipar] = eff_drange1
    #outfile.write("{} {} {} {} {}\n".format(parval,drange1,INL1,avdcnoise1,eff_drange1))
    del tf
    gc.collect()

outfile.close()

#now plot
# fig = plt.figure("Dynamic Range")
# ax = plt.subplot(111)
# ax.set_xlabel(parname, fontsize=14)
# ax.set_ylabel("dynamic range (mV)", fontsize=14)
# ax.plot(parvals,drange,marker='o',linewidth=0,color='k')
# fig.savefig(datadir1+"drange.png")

# fig = plt.figure("DC Noise")
# ax = plt.subplot(111)
# ax.set_xlabel(parname, fontsize=14)
# ax.set_ylabel("DC noise (mV)", fontsize=14)
# ax.set_ylim(0,4)
# ax.plot(parvals,avdcnoise,marker='o',linewidth=0,color='k')
# fig.savefig(datadir1+"dcnoise.png")

# fig = plt.figure("INL")
# ax = plt.subplot(111)
# ax.set_xlabel(parname, fontsize=14)
# ax.set_ylabel("Integrated Non Linearity (ADC counts)", fontsize=14)
# ax.plot(parvals,INL,marker='o',linewidth=0,color='k')
# fig.savefig(datadir1+"inl.png")

# fig = plt.figure("Effective Dynamic Range")
# ax = plt.subplot(111)
# ax.set_xlabel(parname, fontsize=14)
# ax.set_ylabel("effective dynamic range (bit)", fontsize=14)
# ax.plot(parvals,eff_drange,marker='o',linewidth=0,color='k')
# fig.savefig(datadir1+"effdrane.png")
