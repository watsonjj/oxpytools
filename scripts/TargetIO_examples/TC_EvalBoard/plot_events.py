import target_driver
import target_io
import matplotlib.pyplot as plt
import numpy as np

filename = "/home/cta/luigi/TC_EvalBoard/data/20160401/test/Isel_2400/TFdata_VpedDAC_1575.fits"
Nsamples = 96

reader = target_io.EventFileReader(filename)

NEvents  = reader.GetNEvents()

ax = plt.subplot(111)

ampl = np.zeros([NEvents,Nsamples])

for ievt in range(NEvents):
    rawdata = reader.GetEventPacket(ievt,0)#second entry is asic/channel, only one enabled in this case
    packet = target_driver.DataPacket()
    packet.Assign(rawdata, reader.GetPacketSize())
    wf = packet.GetWaveform(0)
    for sample in range(Nsamples):
        ampl[ievt,sample] = wf.GetADC(sample+12)

    ax.plot(ampl[ievt])

ax.set_xlabel("sample",fontsize=14)
ax.set_ylabel("raw ADC counts", fontsize=14)

plt.show()
