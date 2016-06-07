import target_io
# import target_driver
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

filename = "/Users/Jason/Software/outputs/libCHEC/sky/cameradata_run1594.fits"
nSamples = 96
nchannel = 16
nasic = 4

reader = target_io.EventFileReader(filename)

nEvents = reader.GetNEvents()

print("number of events", nEvents)

# ampl = np.zeros([nEvents,nasic, nchannel, nSamples])
# for asic in range(nasic):
# #for asic in range(1):
#         fig = plt.figure('ASIC {}'.format(asic), (20.,8.))
#         fig.suptitle("Asic {}".format(asic), fontsize=16, fontweight='bold', bbox={'boxstyle':'round4','facecolor':'gray', 'alpha':0.5})
#         gs = gridspec.GridSpec(4,4)
#         gs.update(left=0.06,right= 0.98, top=0.93, bottom = 0.05, wspace=0.25,hspace=0.25)
#
#         for channel in range(nchannel):
#
#             ax1 = plt.subplot(gs[int(channel/4), int(channel%4)])
#             ax1.set_xlabel('ns')
#             ax1.set_ylabel('raw ADC counts')
#             ax1.grid(True)
#             plt.text(0.65, 0.9," Channel {}".format(channel), transform = ax1.transAxes, bbox={'boxstyle':'roundtooth','facecolor':'grey', 'alpha':0.2})
#
#             for ievt in range(nEvents):
#                 rawdata = reader.GetEventPacket(ievt,12*64+asic*nchannel+channel)
#                 packet = target_driver.DataPacket()
#                 packet.Assign(rawdata, reader.GetPacketSize())
#                 wf = packet.GetWaveform(0)
#                 # print(dir(wf))
#                 # print(wf.GetADC())
#
#
#                 #print asic, wf.GetASIC()
#                 #print channel, wf.GetChannel()
#                 for sample in range(nSamples):
#                     ampl[ievt,asic, channel, sample] = wf.GetADC(sample)
#                 # if np.any(ampl == 0):
#                 #    print np.where(ampl == 0)
#                 ax1.plot(ampl[0, asic, channel],alpha=1, linewidth = 0.1)
#
#         #fig.savefig(basename+"_ADCwfAsic{}.pdf".format(asic))
#
# plt.show()

event_index = 0
for ipack in range(reader.GetNPacketsPerEvent()):
    event_packet = reader.GetEventPacket(event_index, ipack)
    packet = target_io.DataPacket()
    packet.Assign(event_packet, reader.GetPacketSize())
    if packet.GetSlotID() > 0:
        continue
    print("SlotID = ", packet.GetSlotID())
    print("ChannelID = ", packet.GetChannelID())
    print("ASICID = ", packet.GetASICID())
    # print("WaveformLength = ", packet.GetWaveformLength())
    # print("NumberOfWaveforms = ", packet.GetNumberOfWaveforms())
    # print(dir(packet))
    for iwav in range(packet.GetNumberOfWaveforms()):
        wav = packet.GetWaveform(iwav)
        print("ChannelID = ", wav.GetChannel())
        print("ASIC = ", wav.GetASIC())
        waveform = np.zeros(wav.GetSamples())
        # print("Samples = ", wav.GetSamples())
        for isam in range(wav.GetSamples()):
            waveform[isam] = wav.GetADC(isam)

        plt.plot(waveform)
        # plt.show()

