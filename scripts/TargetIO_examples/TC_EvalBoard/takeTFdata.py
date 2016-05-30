import target_driver
import target_io
import time
import numpy as np
import pdb
import os

outdir = "/home/cta/luigi/TC_EvalBoard/data/20160330/test7/"
os.mkdir(outdir)

Vped_DAC = np.arange(495,4096,1)
#Vped_DAC = np.arange(495,499,2)

my_ip = "192.168.0.1"
board_ip = "192.168.0.123"

board_def = "/home/cta/TARGET/TargetDriver/config/TEC_FPGA_Firmware0xFEDA0002.def"
asic_def = "/home/cta/TARGET/TargetDriver/config/TEC_ASIC.def"

board = target_driver.TargetModule(board_def, asic_def, 0)
board.ConnectToServer(my_ip, 8200, board_ip, 8105)

#select alternative Wilkinson clock freq
board.WriteSetting("WilkinsonClockFreq",0b01)
#make V boundaries of ramp larger
board.WriteASICSetting('SBias',0,1050)
#shorten ramp
board.WriteSetting("RampSignalDuration",510)
#increase slew rate
board.WriteASICSetting("Isel",0,2180)
board.WriteASICSetting("PUbias",0,3100)
board.WriteASICSetting("CMPbias",0,1440)
board.WriteASICSetting("CMPbias2",0,600)

#Enable synchronous trigger mode
board.WriteSetting("TriggerDelay", 100)
board.WriteSetting("TACK_TriggerType", 0x0)
board.WriteSetting("TACK_TriggerMode", 0x0)
board.WriteSetting("TACK_EnableTrigger", 0x10000)
#Enable one channel (Channel 0) and enable zero supression
board.WriteSetting("EnableChannel0", 0x1)
board.WriteSetting("EnableChannel1", 0x0)
board.WriteSetting("Zero_Enable", 0x1)
board.WriteSetting("NumberOfBuffers", 3)
# last block not usable according to Manuel note
kNPacketsPerEvent = 1#reading only one channel
#if multiple channels depend on FPGA settings
kPacketSize = 278#you can read it on wireshark
# check for data size in bytes
kBufferDepth = 1000

for Vped in Vped_DAC:

    board.WriteSetting("Vped_value", int(Vped))
    time.sleep(1)
    board.WriteSetting("ExtTriggerDirection", 0x1)  # switch on sync trigger

    listener = target_io.DataListener(kBufferDepth, kNPacketsPerEvent, kPacketSize)
    target_driver.TargetModule.DataPortPing(my_ip, board_ip)#akira's hack
    listener.AddDAQListener(my_ip)
    listener.StartListening()

    outfile = outdir+"TFdata_VpedDAC_{:04d}.fits".format(Vped)
    print "taking data for Vped {} DAC".format(Vped)
    writer = target_io.EventFileWriter(outfile, kNPacketsPerEvent, kPacketSize)
    buf = listener.GetEventBuffer()
    writer.StartWatchingBuffer(buf)

    time.sleep(2)
    pdb.set_trace()

    writer.StopWatchingBuffer()
    board.WriteSetting("ExtTriggerDirection", 0x0)

    buf.Flush()
    writer.Close()

board.CloseSockets()


# for Vped in Vped_DAC:
#     board.WriteSetting("Vped_value", int(Vped))
#     time.sleep(0.5)
#     outfile = outdir+"TFdata_VpedDAC_{:04d}.fits".format(Vped)
#     print "taking data for Vped {} DAC".format(Vped)
#
#     writer = target_io.EventFileWriter(outfile, kNPacketsPerEvent, kPacketSize)
#     writer.StartWatchingBuffer(buf)
#
#     time.sleep(2)
#
#     writer.StopWatchingBuffer()
#     time.sleep(0.5)
#     writer.Close()






