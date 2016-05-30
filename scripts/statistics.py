from ctapipe.io.hessio import hessio_event_source
from ctapipe.utils.datasets import get_path

# input_filepath = "/Users/Jason/Software/outputs/sim_telarray/simtel_run4_gcts_hnsb.gz"
# input_filepath = "/Users/Jason/Software/outputs/sim_telarray/gamma_20deg_0deg_run10251___cta-prod3-merged_desert-2150m-Paranal-subarray-3_cone10.simtel.gz"
# input_filepath = "/Users/Jason/Software/outputs/sim_telarray/simtel_runmeudon_gamma_30tel_30deg_1.gz"
input_filepath = "/Users/Jason/Software/outputs/sim_telarray/meudon_proton/simtel_runmeudon_proton_30tel_30deg_1.gz"


source = hessio_event_source(get_path(input_filepath), max_events=100)

event = next(source)
tels = list(event.dl0.tels_with_data)

print(tels)

Qt = event.mc.tel[tels[0]].photo_electrons
if Qt.any() > 0:
    print("File has true pe")
else:
    print("No true pe data")