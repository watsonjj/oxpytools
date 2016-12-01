import numpy as np
from targetpipe.io.files import CHECInputFile as InputFile

from IPython import embed


def plot_pulse_shape():
    pass

if __name__ == '__main__':
    input_filepath = "/Users/Jason/Software/outputs/sim_telarray/simtel_run4_gcts_hnsb.gz"
    input_file = InputFile(input_filepath)
    source = input_file.read()

    old_refshapes = [np.empty(1), np.empty(1)]
    old_refstep = None
    old_lrefshape = None
    old_time_slice = None

    event = next(source)
    tels = list(event.dl0.tels_with_data)
    t = tels[0]

    refshapes = event.mc.tel[t].refshapes[0]
    refstep = event.mc.tel[t].refstep
    lrefshape = event.mc.tel[t].lrefshape
    time_slice = event.mc.tel[t].time_slice



    embed()

    # for waveforms in source:
        # tels = list(waveforms.dl0.tels_with_data)
        # for t in tels:
        #     print(t)
        #
        #     refshapes = [np.empty(1), np.empty(1)]
        #
        #     for chan in range(waveforms.dl0.tel[t].num_channels):
        #         refshapes[chan] = waveforms.mc.tel[t].refshapes[chan]
        #         print(np.array_equal(refshapes[chan], old_refshapes[chan]))
        #
        #     refstep = waveforms.mc.tel[t].refstep
        #     lrefshape = waveforms.mc.tel[t].lrefshape
        #     time_slice = waveforms.mc.tel[t].time_slice
        #
        #     print(refstep)
        #     print(lrefshape)
        #     print(time_slice)
        #     print(" ")
        #
        #     old_refshapes = refshapes
        #     old_refstep = refstep
        #     old_lrefshape = lrefshape
        #     old_time_slice = time_slice