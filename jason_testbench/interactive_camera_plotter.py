import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Slider
from matplotlib import pyplot as plt

import tkinter as tk

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
    def __init__(self, waveform_data, geometry,
                 camera_axis, waveform_axis_list):
        self.__waveform_data = None
        self.__waveform_axis_list = None

        self.camera_axis = camera_axis
        self.waveform_axis_list = []
        self.waveform_title_list = []
        self.waveforms = []
        self.pixel_switch_num = 0

        self.colors = itertools.cycle(['r', 'b', 'c', 'm', 'y', 'k', 'w', 'g'])
        self.current_color = 'r'

        self.active_pixels = []
        self.active_pixel_patches = []
        self.active_pixel_labels = []

        CameraDisplay.__init__(self, geometry, ax=camera_axis)
        self._active_pixel.set_linewidth(1.5)

        self.geom = geometry
        self.waveform_data = waveform_data
        self.waveform_axis_list = waveform_axis_list

        geometry_text = "Geometry = {}".format(geometry.cam_id)
        self.camera_axis.text(0.01, 0.99, geometry_text,
                              horizontalalignment='left',
                              verticalalignment='top',
                              transform=self.camera_axis.transAxes)

    @property
    def waveform_data(self):
        return self.__waveform_data

    @waveform_data.setter
    def waveform_data(self, ndarray):
        self.__waveform_data = ndarray
        for waveform_index in range(len(self.waveforms)):
            pixid = self.active_pixels[waveform_index]
            self.update_waveform(waveform_index, pixid)

    @property
    def waveform_axis_list(self):
        return self.__waveform_axis_list

    @waveform_axis_list.setter
    def waveform_axis_list(self, wlist):
        self.__waveform_axis_list = wlist
        for ax in wlist:
            self.current_color = next(self.colors)
            data = self.waveform_data[0, :]
            waveform, = ax.plot(data)
            waveform.set_color(self.current_color)
            ax.set_xlim(0, data.size - 1)
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("Amplitude (ADC)")
            title = ax.set_title("Waveform (pixel: {})".format(0))
            title.set_color(self.current_color)

            self.waveform_title_list.append(title)
            self.waveforms.append(waveform)

            self._add_pixel_annotation()

    def update_waveform(self, waveform_index, new_pix_id):
        waveform_axis = self.waveform_axis_list[waveform_index]
        waveform_title = self.waveform_title_list[waveform_index]
        waveform = self.waveforms[waveform_index]
        data = self.waveform_data[new_pix_id, :]
        waveform.set_xdata(np.arange(data.size))
        waveform.set_ydata(data)
        waveform_axis.relim()  # make sure all the data fits
        waveform_axis.autoscale()
        waveform_axis.set_xlim(0, data.size - 1)
        waveform_title.set_text("Waveform (pixel: {})".format(new_pix_id))
        plt.draw()

    def _add_pixel_annotation(self):
        self.active_pixels.append(0)
        self.active_pixel_patches.append(copy.copy(self._active_pixel))
        self.active_pixel_patches[-1].set_facecolor('none')
        self.active_pixel_patches[-1].set_edgecolor(self.current_color)
        self.camera_axis.add_patch(self.active_pixel_patches[-1])
        new_pixel_label = self.camera_axis.text(0, 0, 0)
        new_pixel_label.update_from(self._active_pixel_label)
        self.active_pixel_labels.append(new_pixel_label)

    def _on_pick(self, event):
        """ handler for when a pixel is clicked """
        pixid = event.ind[-1]
        n = len(self.waveforms)
        pixel_index = self.pixel_switch_num % n if n > 0 else 0
        if not self.active_pixel_patches:
            self._add_pixel_annotation()
        self.set_active_pixel(pixel_index, pixid)

    def set_active_pixel(self, pixel_index, pix_id):
        log.info('[plot] Switching active pixel (pixel={})'.format(pix_id))
        pixel_patch = self.active_pixel_patches[pixel_index]
        pixel_label = self.active_pixel_labels[pixel_index]

        xx, yy, aa = u.Quantity(self.geom.pix_x[pix_id]).value, \
                     u.Quantity(self.geom.pix_y[pix_id]).value, \
                     u.Quantity(np.array(self.geom.pix_area)[pix_id])
        rr = sqrt(aa)
        if self.geom.pix_type.startswith("hex"):
            pixel_patch.xy = (xx, yy)
        else:
            pixel_patch.xy = (xx - rr / 2., yy - rr / 2.)
        pixel_patch.set_visible(True)
        pixel_label.set_x(xx)
        pixel_label.set_y(yy + rr)
        pixel_label.set_text("{:003d}".format(pix_id))
        pixel_label.set_visible(True)

        if self.waveforms:
            n = len(self.waveforms)
            waveform_index = self.pixel_switch_num % n
            self.active_pixels[waveform_index] = pix_id
            self.update_waveform(waveform_index, pix_id)

        self.pixel_switch_num += 1
        self.update()


class InteractiveSourcePlotter:
    def __init__(self, input_file, fig, camera_ax, waveform_ax_list):
        self.__input_file = None
        self.__event_id_list = None
        self.__event = None
        self.__event_index = None
        self.__event_id = None
        self.__telid = None
        self.__geom = None
        self.__waveform_data = None
        self.__camera_image = None
        self.__tels_with_data = None
        self.__current_time = None

        # self.event = None
        # self.event_index = None
        # self.event_id = None
        # self.telid = None

        self.waveform_data_type = 'adc'
        self.camera_image_type = 'adc'
        # self.waveform_data_type = 'cells'
        # self.camera_image_type = 'cells'
        self.data_min = None
        self.data_max = None
        self.camera = None
        self.ts_list = None

        self.geom_dict = None
        self.colorbar_added = False
        self.colorbar = None

        self._root = tk.Tk()
        self._root.wm_title("Embedding in TK")
        self._canvas = FigureCanvasTkAgg(fig, master=self._root)
        self._canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.tk_activated = False
        self.s_event_index = None
        self.s_event_id = None
        self.s_telid = None

        # self.fig = fig
        self.camera_ax = camera_ax
        self.waveform_ax_list = waveform_ax_list
        self.input_file = input_file

    @property
    def input_file(self):
        return self.__input_file

    @input_file.setter
    def input_file(self, file):
        self.__input_file = file
        self.geom_dict = {}
        self.event_id_list = file.get_list_of_event_ids()
        self.event = file.get_event(0)

    @property
    def event_id_list(self):
        return self.__event_id_list

    @event_id_list.setter
    def event_id_list(self, elist):
        self.__event_id_list = elist
        if self.s_event_index:
            self.s_event_index.config(to=len(elist))
        if self.s_event_id:
            self.s_event_index.config(values=tuple(elist))

    @property
    def event(self):
        return self.__event

    @event.setter
    def event(self, event_container):
        self.__event = event_container
        self.__event_index = event_container.count
        self.__event_id = event_container.dl0.event_id
        log.info('[plot] Switching event (index={}, id={})'
                 .format(self.event_index, self.event_id))
        self.tels_with_data = list(event_container.dl0.tels_with_data)
        if self.telid not in self.tels_with_data:
            self.telid = self.tels_with_data[0]
        else:
            # Reload all elements (using same telescope)
            self.telid = self.__telid

    @property
    def event_index(self):
        return self.__event_index

    @event_index.setter
    def event_index(self, val):
        self.event = self.input_file.get_event(val)

    @property
    def event_id(self):
        return self.__event_id

    @event_id.setter
    def event_id(self, val):
        self.event = self.input_file.get_event(val, True)

    @property
    def tels_with_data(self):
        return self.__tels_with_data

    @tels_with_data.setter
    def tels_with_data(self, tlist):
        self.__tels_with_data = tlist
        if self.s_telid:
            self.s_telid.config(values=tuple(tlist))

    @property
    def telid(self):
        return self.__telid

    @telid.setter
    def telid(self, val):

        # If tel id has changed delete camera
        if not self.__telid == val:
            self.camera_ax.cla()
            for ax in self.waveform_ax_list:
                ax.cla()
            self.camera = None
        self.__telid = val

        log.info('[plot] Switching telescope (telid={})'.format(self.telid))

        self.set_waveform_data()

        if not self.camera:
            pixel_pos = self.event.meta.pixel_pos[val]
            optical_foclen = self.event.meta.optical_foclen[val]
            if val not in self.geom_dict:
                self.geom = CameraGeometry.guess(*pixel_pos, optical_foclen)
                self.geom_dict[val] = self.geom
            else:
                self.geom = self.geom_dict[val]

        self.current_time = 0

        self.setup_matplotlib_widgets()

        self._canvas.show()
        if self.tk_activated:
            self.update_tk_elements()

    @property
    def geom(self):
        return self.__geom

    @geom.setter
    def geom(self, geometry):
        self.__geom = geometry

        self.camera = EventViewer(self.waveform_data, self.geom,
                                  self.camera_ax, self.waveform_ax_list)
        self.camera.cmap = plt.cm.viridis
        self.camera.enable_pixel_picker()

        if not self.colorbar:
            self.camera.add_colorbar(ax=self.camera_ax,
                                     label="Amplitude (ADC)")
            self.colorbar = self.camera.colorbar
        else:
            self.camera.colorbar = self.colorbar
            self.camera.update(True)

    @property
    def waveform_data(self):
        return self.__waveform_data

    @waveform_data.setter
    def waveform_data(self, ndarray):
        self.__waveform_data = ndarray
        if self.camera:
            self.camera.waveform_data = ndarray

    @property
    def camera_image(self):
        return self.__camera_image

    @camera_image.setter
    def camera_image(self, array):
        self.__camera_image = array
        if self.camera:
            self.camera.image = array
            self.camera.set_limits_minmax(self.data_min, self.data_max)

    @property
    def current_time(self):
        return self.__current_time

    @current_time.setter
    def current_time(self, val):
        self.__current_time = val
        log.info('[plot] Switching current time (time={:.2f})'
                 .format(self.current_time))
        if self.ts_list:
            for ts in self.ts_list:
                if not val == ts.val:
                    xy = ts.poly.xy
                    xy[2] = val, 1
                    xy[3] = val, 0
                    ts.poly.xy = xy
                    ts.valtext.set_text(ts.valfmt % val)
                    if ts.drawon:
                        ts.ax.figure.canvas.draw_idle()
                    ts.val = val
        self.set_camera_image()

    def set_waveform_data(self):
        if self.waveform_data_type == 'adc':
            array = self.event.dl0.tel[self.telid].adc_samples[0]
            self.data_min = np.min(array[np.nonzero(array)])
            self.data_max = np.max(array)
            self.waveform_data = array
        elif self.waveform_data_type == 'cells':
            array = self.event.dl0.tel[self.telid].cells
            self.waveform_data = array

    def set_camera_image(self):
        if self.camera_image_type == 'adc':
            self.camera_image = self.waveform_data[:, int(self.current_time)]
        elif self.waveform_data_type == 'cells':
            array = self.waveform_data[:, int(self.current_time)]
            self.data_min = np.min(array[np.nonzero(array)])-5
            self.data_max = np.max(array)
            self.camera_image = array

    def setup_matplotlib_widgets(self):
        if not self.ts_list:
            log.info('[plot] Adding matplotlib widgets')
            self.ts_list = []
            for ax in self.waveform_ax_list:
                time_slider_ax = plt.axes(ax.get_position())
                time_slider_ax.patch.set_alpha(0)
                time_slider = Slider(time_slider_ax, '', *ax.get_xlim(),
                                     valfmt='%0.0f')
                time_slider.poly.set_color('r')
                time_slider.poly.set_facecolor('none')
                time_slider.vline.set_visible(False)
                self.ts_list.append(time_slider)

            def update_time_slider(val):
                self.current_time = val

            for ts in self.ts_list:
                ts.on_changed(update_time_slider)
        else:
            # Update lims when tel is changed
            for ts, ax in zip(self.ts_list, self.waveform_ax_list):
                ts.valmin = ax.get_xlim()[0]
                ts.valmax = ax.get_xlim()[1]
                ts.ax.set_xlim(*ax.get_xlim())

    # def add_check_buttons(self):

    def add_tk_elements(self):
        log.info('[plot] Adding tkinter elements')
        self.tk_activated = True
        self.add_tk_event_index_picker()
        self.add_tk_event_id_picker()
        self.add_tk_telid_picker()
        self.add_tk_print_values()

        self.add_quit()

    def add_tk_event_index_picker(self):
        l_event = tk.Label(self._root, text="Event_index")
        l_event.pack(side=tk.LEFT)
        self.s_event_index = tk.Spinbox(self._root, width=6, from_=0,
                                        to=len(self.event_id_list) - 1)
        self.s_event_index.pack(side=tk.LEFT)

        def callback():
            val = round(float(self.s_event_index.get()))
            all_vals = np.arange(len(self.event_id_list))
            event_req = all_vals[np.abs(all_vals - val).argmin()]
            self.event_index = event_req

        b = tk.Button(self._root, text="get", width=3, command=callback)
        b.pack(side=tk.LEFT)

    def add_tk_event_id_picker(self):
        l_event = tk.Label(self._root, text="Event_id")
        l_event.pack(side=tk.LEFT)
        self.s_event_id = tk.Spinbox(self._root, width=6,
                                     values=tuple(self.event_id_list))
        self.s_event_id.pack(side=tk.LEFT)

        def callback():
            val = round(float(self.s_event_id.get()))
            all_vals = np.array(self.event_id_list)
            event_req = all_vals[np.abs(all_vals - val).argmin()]
            self.event_id = event_req

        b = tk.Button(self._root, text="get", width=3, command=callback)
        b.pack(side=tk.LEFT)

    def add_tk_telid_picker(self):
        l_event = tk.Label(self._root, text="telid")
        l_event.pack(side=tk.LEFT)
        self.s_telid = tk.Spinbox(self._root, width=6,
                                  values=tuple(self.tels_with_data))
        self.s_telid.pack(side=tk.LEFT)

        def callback():
            val = round(float(self.s_telid.get()))
            all_vals = np.array(self.tels_with_data)
            tel_req = all_vals[np.abs(all_vals - val).argmin()]
            self.telid = tel_req

        b = tk.Button(self._root, text="get", width=3, command=callback)
        b.pack(side=tk.LEFT)

    def add_tk_print_values(self):
        def callback():
            for pix in self.camera.active_pixels:
                val = self.waveform_data[pix, int(self.current_time)]
                log.info('[plot] Event_index={}, Event_id={}, telid={}, '
                         'Pixel={}, Time={:.2f}, Value={}'
                         .format(self.event_index, self.event_id, self.telid,
                                 pix, self.current_time, val))

        b = tk.Button(self._root, text="Print", width=5, command=callback)
        b.pack(side=tk.LEFT)

    def add_quit(self):
        import sys
        button = tk.Button(master=self._root, text='Quit', command=sys.exit)
        button.pack(side=tk.BOTTOM)

    def update_tk_elements(self):
        self.s_event_index.delete(0, "end")
        self.s_event_index.insert(0, self.event_index)
        self.s_event_id.delete(0, "end")
        self.s_event_id.insert(0, self.event_id)
        self.s_telid.delete(0, "end")
        self.s_telid.insert(0, self.telid)


def main():
    fp = "/Users/Jason/Software/outputs/libCHEC/sky/cameradata_run1594.fits"
    reader = 'targetio'
    # fp = "/Users/Jason/Software/ctapipe/ctapipe-extra/datasets/gamma_test.simtel.gz"
    # reader = 'hessio'

    input_file = InputFile(fp, reader)

    fig = plt.figure(figsize=(15, 7))
    ax0 = fig.add_subplot(1, 2, 1)
    ax1 = fig.add_subplot(2, 2, 2)
    ax2 = fig.add_subplot(2, 2, 4)
    fig.subplots_adjust(hspace=.5)

    plot = InteractiveSourcePlotter(input_file, fig, ax0, [ax1, ax2])

    plot.add_tk_elements()
    tk.mainloop()

if __name__ == '__main__':
    main()
