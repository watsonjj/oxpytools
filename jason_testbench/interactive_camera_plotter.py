import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Slider
from matplotlib import pyplot as plt


import tkinter as Tk

import itertools
from ctapipe.io import CameraGeometry
from ctapipe.visualization import CameraDisplay
from targetpipe.io.files import CHECInputFile as InputFile

from astropy import units as u
import numpy as np
from numpy import sqrt
import copy
from astropy import log



class EventViewer(CameraDisplay):
    def __init__(self, event, camera_axis):
        self.event = event
        self.event_tels = list(event.dl0.tels_with_data)
        self.camera_axis = camera_axis
        self.waveform_axis_list = []
        self.waveform_title_list = []
        self.waveforms = []
        self.pixel_switch_num = 0

        self.telid = 0
        self.geom = CameraGeometry.guess(*event.meta.pixel_pos[self.telid],
                                         event.meta.optical_foclen[self.telid])
        self.data = event.dl0.tel[self.telid].adc_samples[0]

        CameraDisplay.__init__(self, self.geom, ax=camera_axis)

        self.colors = itertools.cycle(['r', 'b', 'c', 'm', 'y', 'k', 'w', 'g'])
        self.current_color = 'r'

        self.active_pixels = []
        self.active_pixel_labels = []

    def add_pixel_annotation(self):
        self.active_pixels.append(copy.copy(self._active_pixel))
        self.active_pixels[-1].set_facecolor(self.current_color)
        self.camera_axis.add_patch(self.active_pixels[-1])
        new_pixel_label = self.camera_axis.text(0, 0, 0)
        new_pixel_label.update_from(self._active_pixel_label)
        self.active_pixel_labels.append(new_pixel_label)

    def add_waveform_display(self, ax=None):
        self.current_color = next(self.colors)

        axes = ax if ax is not None else plt.gca()
        waveform, = axes.plot(self.data[0, :])
        waveform.set_color(self.current_color)
        axes.set_xlabel("Time (ns)")
        axes.set_ylabel("Amplitude (ADC)")
        title = axes.set_title("Waveform (pixel: {})".format(0))
        title.set_color(self.current_color)

        self.waveform_axis_list.append(axes)
        self.waveform_title_list.append(title)
        self.waveforms.append(waveform)

        self.add_pixel_annotation()

        return waveform

    def switch_waveform(self, pix_id):
        n = len(self.waveforms)
        waveform_axis = self.waveform_axis_list[self.pixel_switch_num % n]
        waveform_title = self.waveform_title_list[self.pixel_switch_num % n]
        waveform = self.waveforms[self.pixel_switch_num % len(self.waveforms)]
        waveform.set_ydata(self.data[pix_id, :])
        waveform_axis.relim()  # make sure all the data fits
        waveform_axis.autoscale()
        waveform_title.set_text("Waveform (pixel: {})".format(pix_id))
        plt.draw()

    def _on_pick(self, event):
        """ handler for when a pixel is clicked """
        n = len(self.waveforms)
        i = self.pixel_switch_num % n if n > 0 else 0
        if not self.active_pixels:
            self.add_pixel_annotation()
        pixel = self.active_pixels[i]
        pixel_label = self.active_pixel_labels[i]

        pix_id = event.ind[-1]
        xx, yy, aa = u.Quantity(self.geom.pix_x[pix_id]).value, \
                     u.Quantity(self.geom.pix_y[pix_id]).value, \
                     u.Quantity(np.array(self.geom.pix_area)[pix_id])
        if self.geom.pix_type.startswith("hex"):
            pixel.xy = (xx, yy)
        else:
            rr = sqrt(aa)
            pixel.xy = (xx - rr / 2., yy - rr / 2.)
        pixel.set_visible(True)
        pixel_label.set_x(xx)
        pixel_label.set_y(yy)
        pixel_label.set_text("{:003d}".format(pix_id))
        pixel_label.set_visible(True)
        self.on_pixel_clicked(pix_id)  # call user-function
        self.pixel_switch_num += 1
        self.update()

    def on_pixel_clicked(self, pix_id):
        if self.waveforms:
            self.switch_waveform(pix_id)


class InteractiveSourcePlotter:
    def __init__(self, input_file, canvas, camera_ax, waveform_ax_list):
        self.__input_file = None
        self.__event = None
        self.__event_index = None
        self.__event_id = None
        self.__telid = None

        # self.event = None
        # self.event_index = None
        # self.event_id = None
        # self.telid = None

        self.tels_with_data = None
        self.data = None
        self.data_min = None
        self.data_max = None
        self.camera = None
        self.ts_list = None

        self.geom_dict = None
        self.colorbar_added = False
        self.colorbar = None

        self.canvas = canvas
        self.camera_ax = camera_ax
        self.waveform_ax_list = waveform_ax_list
        self.input_file = input_file

    @property
    def input_file(self):
        return self.__input_file

    @input_file.setter
    def input_file(self, file):
        self.__input_file = file
        self.geom_dict = None
        self.event = file.get_event(0)

    @property
    def event(self):
        return self.__event

    @event.setter
    def event(self, event_container):
        print("event")
        self.__event = event_container
        self.__event_index = event_container.count
        self.__event_id = event_container.dl0.event_id
        self.tels_with_data = list(event_container.dl0.tels_with_data)
        self.telid = self.tels_with_data[0]

    @property
    def event_index(self):
        return self.__event_index

    @event_index.setter
    def event_index(self, val):
        print("event_index")
        self.event = self.input_file.get_event(val)

    @property
    def event_id(self):
        return self.__event_id

    @event_id.setter
    def event_id(self, val):
        print("event_id")
        self.event = self.input_file.get_event(val, True)

    @property
    def telid(self):
        return self.__telid

    @telid.setter
    def telid(self, val):
        self.__telid = val

        self.camera_ax.cla()
        for ax in self.waveform_ax_list:
            ax.cla()

        self.data = self.event.dl0.tel[val].adc_samples[0]
        self.data_min = np.min(self.data[np.nonzero(self.data)])
        self.data_max = np.max(self.data)

        self.camera = EventViewer(self.event, self.camera_ax)
        for ax in self.waveform_ax_list:
            self.camera.add_waveform_display(ax)
        self.camera.image = self.data[:, 0]
        self.camera.cmap = plt.cm.viridis
        self.camera.set_limits_minmax(self.data_min * 1.05,
                                      self.data_max * 0.95)
        self.camera.enable_pixel_picker()

        print('colorbar', self.colorbar)
        if not self.colorbar:
            self.camera.add_colorbar(ax=self.camera_ax, label="Amplitude (ADC)")
            self.colorbar = self.camera.colorbar
        else:
            self.camera.colorbar = self.colorbar
            self.camera.update(True)

        if not self.ts_list:
            self.ts_list = self.add_time_sliders()

        self.canvas.show()

    def add_time_sliders(self):
        ts_list = []
        for ax in self.waveform_ax_list:
            time_slider_ax = plt.axes(ax.get_position())
            time_slider_ax.patch.set_alpha(0)
            time_slider = Slider(time_slider_ax, '', *ax.get_xlim(), valfmt='%0.0f')
            time_slider.poly.set_color('r')
            time_slider.poly.set_facecolor('none')
            time_slider.vline.set_visible(False)
            ts_list.append(time_slider)

        def update_time_slider(val):
            max_time = np.size(self.data[0])-1
            if val > max_time:
                val = max_time
            time = int(val)
            self.camera.image = self.data[:, time]
            for ts in ts_list:
                if not val == ts.val:
                    ts.set_val(val)
            plt.draw()

        for ts in ts_list:
            ts.on_changed(update_time_slider)

        return ts_list


def main():
    fp = "/Users/Jason/Software/outputs/libCHEC/sky/cameradata_run1594.fits"
    reader = 'targetio'
    # fp = "/Users/Jason/Software/outputs/sim_telarray/simtel_run4_gcts_hnsb.gz"
    # reader = 'hessio'
    telid = 0
    event_index = 9

    input_file = InputFile(fp, reader)
    event_id_list = input_file.get_list_of_event_ids()

    fig = plt.figure(figsize=(15, 7))
    ax0 = fig.add_subplot(1, 2, 1)
    ax1 = fig.add_subplot(2, 2, 2)
    ax2 = fig.add_subplot(2, 2, 4)
    fig.subplots_adjust(hspace=.5)

    root = Tk.Tk()
    root.wm_title("Embedding in TK")

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

    # event = input_file.get_event(event_index)

    plot = InteractiveSourcePlotter(input_file, canvas, ax0, [ax1, ax2])


    # data = event.dl0.tel[telid].adc_samples[0]
    # data_min = np.min(data[np.nonzero(data)])
    # data_max = np.max(data)
    #
    # camera = EventViewer(event, ax0)
    # camera.add_waveform_display(ax1)
    # camera.add_waveform_display(ax2)
    # camera.image = data[:, 0]
    # camera.cmap = plt.cm.viridis
    # camera.set_limits_minmax(data_min*1.05, data_max*0.95)
    # camera.add_colorbar(ax=ax0, label="Amplitude (ADC)")
    # camera.enable_pixel_picker()
    #
    # ts1, ts2 = add_time_sliders(data, camera, ax1, ax2)
    #
    # canvas.show()
    # plt.show()


    # canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

    # Event picker
    l_event = Tk.Label(root, text="Event")
    l_event.pack(side=Tk.LEFT)
    s_event = Tk.Spinbox(root, from_=0, to=len(event_id_list)-1)
    s_event.pack(side=Tk.LEFT)

    def callback():
        val = round(float(s_event.get()))
        all = np.arange(len(event_id_list))
        event_req = all[np.abs(all-val).argmin()]
        s_event.delete(0, "end")
        s_event.insert(0, event_req)
        plot.event_index = event_req

    b = Tk.Button(root, text="get", width=10, command=callback)
    b.pack(side=Tk.LEFT)


    import sys
    button = Tk.Button(master=root, text='Quit', command=sys.exit)
    button.pack(side=Tk.BOTTOM)

    Tk.mainloop()


if __name__ == '__main__':
    main()