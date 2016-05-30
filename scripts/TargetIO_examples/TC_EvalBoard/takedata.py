import target_driver
import target_io
import time
import pdb

outdir = "/home/cta/luigi/TC_EvalBoard/data/20160330/"
outfile = outdir + "test_dump2.fits"
my_ip = "192.168.0.1"
board_ip = "192.168.0.123"

board_def = "/home/cta/TARGET/TargetDriver/config/TEC_FPGA_Firmware0xFEDA0003.def"
asic_def = "/home/cta/TARGET/TargetDriver/config/TEC_ASIC.def"

board = target_driver.TargetModule(board_def, asic_def, 0)
board.ConnectToServer(my_ip, 8200, board_ip, 8105)

#select alternative Wilkinson clock freq
board.WriteSetting("WilkinsonClockFreq",0b01)
#make V boundaries of ramp larger
board.WriteASICSetting('SBias',0,1020)
#shorten ramp
board.WriteSetting("RampSignalDuration",1080)
#increase slew rate
board.WriteASICSetting("Isel",0,2440)


#Enable synchronous trigger mode
board.WriteSetting("TriggerDelay", 100)
board.WriteSetting("TACK_TriggerType", 0x0)
board.WriteSetting("TACK_TriggerMode", 0x0)
board.WriteSetting("TACK_EnableTrigger", 0x10000)
# #Set Vped
board.WriteSetting("Vped_value", 2000)
#Enable one channel (Channel 0) and enable zero supression
board.WriteSetting("EnableChannel0", 0x1)
board.WriteSetting("EnableChannel1", 0x0)
board.WriteSetting("Zero_Enable", 0x1)


board.WriteSetting("NumberOfBuffers", 3)
# last block not usable according to Manuel note

board.WriteSetting("ExtTriggerDirection", 0x1)  # switch on sync trigger

kNPacketsPerEvent = 1#reading only one channel
#if multiple channels depend on FPGA settings
kPacketSize = 278#you can read it on wireshark
# check for data size in bytes
kBufferDepth = 1000

listener = target_io.DataListener(kBufferDepth, kNPacketsPerEvent, kPacketSize)
target_driver.TargetModule.DataPortPing(my_ip, board_ip)#akira's hack
listener.AddDAQListener(my_ip)
listener.StartListening()

writer = target_io.EventFileWriter(outfile, kNPacketsPerEvent, kPacketSize)
buf = listener.GetEventBuffer()
writer.StartWatchingBuffer(buf)

#time.sleep(2)
pdb.set_trace()

writer.StopWatchingBuffer()
board.WriteSetting("ExtTriggerDirection", 0x0)

board.CloseSockets()
buf.Flush()
writer.Close()



