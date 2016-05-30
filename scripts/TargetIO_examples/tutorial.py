import time

import target_driver
import target_io

kNWavesPerPacket = 8
kNSamplesPerWave = 128

def SendDataPacket(simulator, packet, eventID, packetID):
    tack = eventID
    slotID = 0
    eventSequenceNumber = 0
    quad = 0
    row = 0
    col = 0

    tmp = packetID * kNWavesPerPacket
    detectorID = tmp / 64
    tmp %= 64
    asic = tmp / 16
    tmp %= 16;
    ch_offset = (tmp / kNWavesPerPacket) * kNWavesPerPacket

    packet.FillHeader(kNWavesPerPacket, kNSamplesPerWave, slotID, detectorID,
                      eventSequenceNumber, tack, quad, row, col);

    for i in range(kNWavesPerPacket):
        waveform = packet.GetWaveform(i)
        waveform.SetHeader(asic, ch_offset + i, kNSamplesPerWave, False)

    pid = 0
    ret, pid = packet.GetPacketID()

    simulator.SendDataPacket(packet.GetData(), packet.GetPacketSize())


packet = target_driver.DataPacket(kNWavesPerPacket, kNSamplesPerWave)
    
simulator = target_driver.ModuleSimulator("0.0.0.0")
simulator.Start()

module = target_driver.TargetModule(
    "/Users/oxon/Documents/workspace/TargetDriver/config/TM5_FPGA_Firmware0xFEDA003C.def",
    "/Users/oxon/Documents/workspace/TargetDriver/config/TM5_ASIC.def", 0)
module.EstablishSlowControlLink("0.0.0.0", "0.0.0.0")

kNPacketsPerEvent = 256
kPacketSize = 2084
kBufferDepth = 1000

listener = target_io.DataListener(kBufferDepth, kNPacketsPerEvent, kPacketSize)
listener.AddDAQListener("0.0.0.0")
listener.StartListening()

writer = target_io.EventFileWriter("testEventFile.fits", kNPacketsPerEvent, kPacketSize)

buf = listener.GetEventBuffer()
writer.StartWatchingBuffer(buf)

for eventID in range(100):
    for packetID in range(kNPacketsPerEvent):
        SendDataPacket(simulator, packet, eventID, packetID)

buf.Flush()
