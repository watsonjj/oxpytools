import target_driver
import time
import target_io
import run_control
import BrickPowerSupply
import numpy
#import argparse
import random

#minVped = 125
#maxVped = 4061+1
#stepVped = 164
#VpedList = numpy.concatenate([range(minVped,maxVped,stepVped),range(minVped,maxVped,stepVped),range(minVped,maxVped,stepVped),range(minVped,maxVped,stepVped),range(minVped,maxVped,stepVped)])
#VpedList = range(minVped,maxVped,stepVped)
#VpedList = [2000,2000,2000,2000,2000]
VpedList = [2000]*25

#runDuration = 1
#runDurationList = [1,3,10]

runDurationList = [1,3,10]


remote = 0
Ext = 1
PCifFail = 1
InitOnce = 1

#for Vped in VpedList:
if InitOnce:
	if remote:
		InitFail = 1
		while InitFail != 0:
			bps = BrickPowerSupply.BrickPowerSupply('2200','/dev/tty.usbserial')
			time.sleep(0.5)
			bps.open()
			bps.outputOn()
			#bps.setVoltage(12)
			#bps.setMaxCurrent(2.2)
			time.sleep(30)
		        my_ip = "192.168.1.2"
        		tb_ip = "192.168.1.173"
        		board_ip = "192.168.0.173"
        		board_def = "/Users/pkarn/TargetDriver/config/TM7_FPGA_Firmware0xB0000100.def"
        		asic_def = "/Users/pkarn/TargetDriver/config/TM7_ASIC.def"
        		#hostname = run_control.getHostName()
        		#outdirname = run_control.getDataDirname(hostname)
        		#runID = run_control.incrementRunNumber(outdirname)
        		#outFile = "%srun%d.fits" % (outdirname,runID)
        		#print "Writing to: %s" % outFile
        		board = target_driver.TargetModule(board_def, asic_def, 1)
        		InitFail = board.EstablishSlowControlLink(my_ip, board_ip)
			print(InitFail)
			if InitFail:
				bps.outputOff()
				time.sleep(15)

	elif PCifFail:
		InitFail = 1
                while InitFail != 0:
                        bps = BrickPowerSupply.BrickPowerSupply('2200','/dev/tty.usbserial')
                        time.sleep(0.5)
                        bps.open()
                        #bps.outputOn()
                        #bps.setVoltage(12)
                        #bps.setMaxCurrent(2.2)
                        #time.sleep(30)
                        my_ip = "192.168.1.2"
                        tb_ip = "192.168.1.173"
                        board_ip = "192.168.0.173"
                        board_def = "/Users/pkarn/TargetDriver/config/TM7_FPGA_Firmware0xB0000100.def"
                        asic_def = "/Users/pkarn/TargetDriver/config/TM7_ASIC.def"
                        #hostname = run_control.getHostName()
                        #outdirname = run_control.getDataDirname(hostname)
                        #runID = run_control.incrementRunNumber(outdirname)
                        #outFile = "%srun%d.fits" % (outdirname,runID)
                        #print "Writing to: %s" % outFile
                        board = target_driver.TargetModule(board_def, asic_def, 1)
                        InitFail = board.EstablishSlowControlLink(my_ip, board_ip)
                        print(InitFail)
                        if InitFail:
                                bps.outputOff()
                                time.sleep(15)
				bps.outputOn()
				time.sleep(30)


	else:
		my_ip = "192.168.1.2"
		tb_ip = "192.168.1.173"
		board_ip = "192.168.0.173"
		board_def = "/Users/pkarn/TargetDriver/config/TM7_FPGA_Firmware0xB0000100.def"
		asic_def = "/Users/pkarn/TargetDriver/config/TM7_ASIC.def"
		#hostname = run_control.getHostName()
		#outdirname = run_control.getDataDirname(hostname)
		#runID = run_control.incrementRunNumber(outdirname)
		#outFile = "%srun%d.fits" % (outdirname,runID)
		#print "Writing to: %s" % outFile
		board = target_driver.TargetModule(board_def, asic_def, 1)
		print board.EstablishSlowControlLink(my_ip, board_ip)

	time.sleep(0.1)
	print board.Initialise()

#	hostname = run_control.getHostName()
#	outdirname = run_control.getDataDirname(hostname)
#	runID = run_control.incrementRunNumber(outdirname)
#	outFile = "%srun%d.fits" % (outdirname,runID)
#	print "Writing to: %s" % outFile

	time.sleep(0.1)
	tester = target_driver.TesterBoard()
	tester.Init(my_ip, tb_ip)
	time.sleep(0.1)
	tester.EnableReset(True)
	time.sleep(0.1)
	tester.EnableReset(False)
	time.sleep(0.1)

	##board.WriteSetting("TACK_EnableTrigger",0b010000000000000000)
	##board.WriteSetting("ExtTriggerDirection",0)
	time.sleep(0.5)

for runDuration in runDurationList:

	for Vped in VpedList:

        	hostname = run_control.getHostName()
        	outdirname = run_control.getDataDirname(hostname)
        	runID = run_control.incrementRunNumber(outdirname)
        	outFile = "%srun%d.fits" % (outdirname,runID)
        	print "Writing to: %s" % outFile

		#Set Vped, enable all channels for readout
		board.WriteSetting("Vped_value", Vped)
		board.WriteSetting("EnableChannel0", 0xffffffff)
		board.WriteSetting("EnableChannel1", 0xffffffff)
		board.WriteSetting("Zero_Enable", 0x1)

		board.WriteSetting("NumberOfBuffers", 1)
		board.WriteSetting("TriggerDelay", 8)

		#HV Control
		#Start MAX11616, HV monitoring
		#print board.WriteRegister(0x1F, 0x80000000)
		#print board.WriteSetting("HV_Enable", 0x0)
		#time.sleep(0.1)
		##print board.WriteSetting("HV_Enable", 0x1)
		time.sleep(0.1)
		#Set reference voltage to 4.096V, need to do twice?
		board.WriteRegister(0x20, 0x7300080F)
		board.WriteRegister(0x20, 0x7300080F)
		time.sleep(1)
		#Set low side voltage for all trigger channels to 2.048V
		board.WriteRegister(0x20, 0x8280080F)

		kNPacketsPerEvent = 64 #One per channel
		#if multiple channels depend on FPGA settings
		kPacketSize = 150#you can read it on wireshark
		# check for data size in bytes
		kBufferDepth = 3000

		#Set up IO
		listener = target_io.DataListener(kBufferDepth, kNPacketsPerEvent, kPacketSize)
		#target_driver.TargetModule.DataPortPing(my_ip, board_ip)#akira's hack
		listener.AddDAQListener(my_ip)
		listener.StartListening()

		writer = target_io.EventFileWriter(outFile, kNPacketsPerEvent, kPacketSize)
		buf = listener.GetEventBuffer()

		#time.sleep(random.randint(1,20))
		#time.sleep(5)

		#board.WriteSetting("TACK_EnableTrigger",0x10000)
	        #board.WriteSetting("ExtTriggerDirection",1)


		if Ext:
			tester.EnableExternalTrigger(True)
			writer.StartWatchingBuffer(buf)
			time.sleep(runDuration)
			writer.StopWatchingBuffer()
			tester.EnableExternalTrigger(False)

		else:
			#Enable external trigger
			#tester.EnableExternalTrigger(True)

			tester.EnableSoftwareTrigger(True)

			#board.WriteSetting("Row",5)
			#board.WriteSetting("Column",0)
	
			writer.StartWatchingBuffer(buf)

			#Taking data now

			for i in range(runDuration*300):

				tester.SendSoftwareTrigger()
				time.sleep(1/300.)
	
			#time.sleep(1)

			writer.StopWatchingBuffer()

			#Stop trigger
			#tester.EnableExternalTrigger(False)

			tester.EnableSoftwareTrigger(False)

			##board.WriteSetting("HV_Enable", 0x0)
			#Close everything

		#listener.StopListening()

		#board.CloseSockets()
		#tester.CloseSockets()
		buf.Flush()
		writer.Close()
		if remote:
			bps.outputOff()
		if numpy.size(VpedList) == 1:
			time.sleep(0.1)
		else:
			if remote:
				time.sleep(20)
			else:
				time.sleep(5)
	
board.CloseSockets()
tester.CloseSockets()

