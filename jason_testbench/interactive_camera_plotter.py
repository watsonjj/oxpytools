import matplotlib
from ctapipe.calib.camera.calibrators import calibration_parameters, \
    calibrate_event
from ctapipe.calib.camera.integrators import integrator_dict

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
                 camera_axis, waveform_axis_list, waveform_axis_ylabel):
        self.__waveform_data = None
        self.__waveform_axis_list = None
        self.__integration_window = None

        self.camera_axis = camera_axis
        self.waveform_axis_list = []
        self.waveform_title_list = []
        self.waveform_window_list = []
        self.waveforms = []
        self.pixel_switch_num = 0

        self.colors = itertools.cycle(['r', 'b', 'c', 'm', 'y', 'k', 'w', 'g'])
        self.current_color = 'r'

        self.active_pixels = []
        self.active_pixel_patches = []
        self.active_pixel_labels = []

        self.window_start = None
        self.window_end = None


        CameraDisplay.__init__(self, geometry, ax=camera_axis)
        self._active_pixel.set_linewidth(1.5)

        self.geom = geometry
        self.waveform_yaxis_label = waveform_axis_ylabel
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
            ax.set_ylabel(self.waveform_yaxis_label)
            title = ax.set_title("Waveform (pixel: {})".format(0))
            title.set_color(self.current_color)

            # Integration window lines
            line1, = ax.plot([0, 0], ax.get_ylim(), color='b', visible=False)
            line2, = ax.plot([0, 0], ax.get_ylim(), color='b', visible=False)

            self.waveform_title_list.append(title)
            self.waveforms.append(waveform)
            self.waveform_window_list.append((line1, line2))

            self._add_pixel_annotation()

    @property
    def integration_window(self):
        return self.__integration_window

    @integration_window.setter
    def integration_window(self, ndarray):
        self.__integration_window = ndarray
        length = np.sum(ndarray, axis=1)
        self.window_start = np.argmax(ndarray, axis=1)
        self.window_end = self.window_start + length - 1
        self.update_window_position()

    def update_waveform(self, waveform_index, new_pix_id):
        waveform_axis = self.waveform_axis_list[waveform_index]
        waveform_title = self.waveform_title_list[waveform_index]
        line = self.waveform_window_list[waveform_index]
        waveform = self.waveforms[waveform_index]
        data = self.waveform_data[new_pix_id, :]
        waveform.set_xdata(np.arange(data.size))
        waveform.set_ydata(data)
        waveform_axis.relim()  # make sure all the data fits
        waveform_axis.autoscale()
        waveform_axis.set_xlim(0, data.size - 1)
        waveform_axis.set_ylabel(self.waveform_yaxis_label)
        if self.integration_window is not None:
            line[0].set_xdata((self.window_start[new_pix_id],
                               self.window_start[new_pix_id]))
            line[1].set_xdata((self.window_end[new_pix_id],
                               self.window_end[new_pix_id]))
            line[0].set_ydata(waveform_axis.get_ylim())
            line[1].set_ydata(waveform_axis.get_ylim())
        waveform_title.set_text("Waveform (pixel: {})".format(new_pix_id))
        waveform_axis.figure.canvas.draw()

    def update_window_position(self):
        for lines, pix in zip(self.waveform_window_list, self.active_pixels):
            lines[0].set_xdata((self.window_start[pix],
                                self.window_start[pix]))
            lines[1].set_xdata((self.window_end[pix],
                                self.window_end[pix]))

    def update_window_visibility(self, visible):
        for lines, pix in zip(self.waveform_window_list, self.active_pixels):
            lines[0].set_visible(visible)
            lines[1].set_visible(visible)

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
        self.__calibrate = None
        self.__waveform_view = None
        self.__camera_view = None

        # self.waveforms = None
        # self.event_index = None
        # self.event_id = None
        # self.telid = None

        self.data_min = None
        self.data_max = None
        self.camera = None
        self.ts_list = None
        self.calib_params = None
        self.waveform_views_dict = None
        self.camera_views_dict = None

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
        self.om_waveform_view = None
        self.om_camera_view = None
        self.calib_popout = None
        self.wf_view_var = None
        self.c_view_var = None

        self.waveform_yaxis_label = 'Amplitude (ADC)'
        self.camera_axis_label = 'Amplitude (ADC)'

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
        self.set_calibration_params()
        self.update_views_dict()
        self.event_id_list = file.event_id_list
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
        ev = event_container
        log.info('[plot] Switching waveforms (index={}, id={})'
                 .format(self.event_index, self.event_id))
        if self.calibrate and self.calib_params:
            log.info('[plot] Calibrating waveforms')
            ev = calibrate_event(ev, self.calib_params, self.geom_dict)
        self.__event = ev
        self.__event_index = ev.count
        self.__event_id = ev.dl0.event_id
        self.tels_with_data = list(ev.dl0.tels_with_data)
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

        self.waveform_view = self.waveform_view

        if not self.camera:
            pixel_pos = self.event.inst.pixel_pos[val]
            optical_foclen = self.event.inst.optical_foclen[val]
            if val not in self.geom_dict:
                self.geom = CameraGeometry.guess(*pixel_pos, optical_foclen)
                self.geom_dict[val] = self.geom
            else:
                self.geom = self.geom_dict[val]

        if self.calibrate and self.calib_params:
            windows = self.event.dl1.tel[val].integration_window[0]
            self.camera.integration_window = windows
            self.camera.update_window_visibility(True)
        else:
            self.camera.update_window_visibility(False)

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
                                  self.camera_ax, self.waveform_ax_list,
                                  self.waveform_yaxis_label)
        self.camera.cmap = plt.cm.viridis
        self.camera.enable_pixel_picker()

        if not self.colorbar:
            self.camera.add_colorbar(ax=self.camera_ax,
                                     label=self.camera_axis_label)
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
            self.camera.waveform_yaxis_label = self.waveform_yaxis_label
            self.camera.waveform_data = ndarray

    @property
    def camera_image(self):
        return self.__camera_image

    @camera_image.setter
    def camera_image(self, array):
        self.__camera_image = array
        if self.camera:
            self.camera.colorbar.set_label(self.camera_axis_label)
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
        self.camera_view = self.camera_view

    @property
    def calibrate(self):
        return self.__calibrate

    @calibrate.setter
    def calibrate(self, flag):
        self.__calibrate = flag
        # Re-set waveforms so it can be calibrated
        if flag:
            log.info('[plot] Calibration activated')
            self.event = self.event
        else:
            log.info('[plot] Calibration de-activated')
        self.update_views_dict()

    def set_calibration_params(self, cmdline=None):
        if not cmdline:
            cmdline = []
        try:
            params, unk = calibration_parameters(cmdline,
                                                 self.input_file.origin)
        except KeyError:
            params = None
        self.calib_params = params

    @property
    def waveform_view(self):
        if self.__waveform_view is None:
            return 'adc'
        return self.__waveform_view

    @waveform_view.setter
    def waveform_view(self, view):
        if view not in self.waveform_views_dict:
            view = 'adc'
            self.wf_view_var.set(view)
        if not view == self.__waveform_view:
            log.info('[plot] Changing waveform view ({})'.format(view))
        self.__waveform_view = view
        self.waveform_views_dict[view]()

    @property
    def camera_view(self):
        if self.__camera_view is None:
            return 'adc'
        return self.__camera_view

    @camera_view.setter
    def camera_view(self, view):
        if view not in self.camera_views_dict:
            view = 'adc'
            self.c_view_var.set(view)
        if not view == self.__camera_view:
            log.info('[plot] Changing camera view ({})'.format(view))
        self.__camera_view = view
        self.camera_views_dict[view]()

    def update_views_dict(self):
        self.update_waveform_views_dict()
        self.update_camera_views_dict()

    def update_waveform_views_dict(self):
        # Default views
        def view_adc():
            self.waveform_yaxis_label = 'Amplitude (ADC)'
            array = self.event.dl0.tel[self.telid].adc_samples[0]
            # self.data_min = np.min(array[np.nonzero(array)])
            # self.data_max = np.max(array)
            self.waveform_data = array

        views = {
            'adc': lambda: view_adc()
        }

        # Origin-specific views
        if self.input_file:
            if self.input_file.origin == 'targetio':
                def view_cells():
                    self.waveform_yaxis_label = 'CellID'
                    array = self.event.dl0.tel[self.telid].cells
                    self.waveform_data = array

                views['cells'] = lambda: view_cells()

        # Calibration views
        if self.calibrate and self.calib_params:
            def view_pededestal_subtracted():
                self.waveform_yaxis_label = 'Amplitude (ADC-ped)'
                ev = self.event.dl1.tel[self.telid]
                array = ev.pedestal_subtracted_adc[0]
                self.waveform_data = array

            views['pedestal_subtracted'] = lambda: view_pededestal_subtracted()

        self.waveform_views_dict = views
        if self.tk_activated:
            self.update_tk_om(self.om_waveform_view,
                              self._wf_views_callback,
                              self.waveform_views_dict.keys())

    def update_camera_views_dict(self):
        # Default views
        def view_adc():
            self.camera_axis_label = 'Amplitude (ADC)'
            array = self.event.dl0.tel[self.telid].adc_samples[0]
            self.data_min = np.min(array[np.nonzero(array)])
            self.data_max = np.max(array)
            self.camera_image = array[:, int(self.current_time)]

        views = {
            'adc': lambda: view_adc()
        }

        # Origin-specific views
        if self.input_file:
            if self.input_file.origin == 'targetio':
                def view_cells():
                    self.camera_axis_label = 'CellID'
                    wf = self.event.dl0.tel[self.telid].cells
                    array = wf[:, int(self.current_time)]
                    self.data_min = np.min(array[np.nonzero(array)]) - 5
                    self.data_max = np.max(array)
                    self.camera_image = array

                views['cells'] = lambda: view_cells()

        # Calibration views
        if self.calibrate and self.calib_params:
            def view_pededestal_subtracted():
                self.camera_axis_label = 'Amplitude (ADC-ped)'
                ev = self.event.dl1.tel[self.telid]
                array = ev.pedestal_subtracted_adc[0]
                self.data_min = np.min(array[np.nonzero(array)])
                self.data_max = np.max(array)
                self.camera_image = array[:, int(self.current_time)]

            def view_pe_charge():
                self.camera_axis_label = 'Amplitude (pe)'
                ev = self.event.dl1.tel[self.telid]
                array = ev.pe_charge
                self.data_min = np.min(array[np.nonzero(array)])
                self.data_max = np.max(array)
                self.camera_image = array[:]

            def view_peak_pos():
                self.camera_axis_label = 'Peakpos (ns)'
                ev = self.event.dl1.tel[self.telid]
                array = ev.peakpos
                self.data_min = np.min(array[np.nonzero(array)])
                self.data_max = np.max(array)
                self.camera_image = array[:]

            views['pedestal_subtracted'] = lambda: view_pededestal_subtracted()
            views['pe_charge'] = lambda: view_pe_charge()
            views['peak_pos'] = lambda: view_peak_pos()

        self.camera_views_dict = views
        if self.tk_activated:
            self.update_tk_om(self.om_camera_view,
                              self._c_views_callback,
                              self.camera_views_dict.keys())

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
        self.add_tk_waveform_views_dropdown()
        self.add_tk_camera_views_dropdown()
        self.add_tk_calibration_settings_button()
        self.add_tk_calibrate_checkbox()

        # self.add_quit()

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

    def _wf_views_callback(self, val):
        self.wf_view_var.set(val)
        self.waveform_view = val

    def add_tk_waveform_views_dropdown(self):
        choices = self.waveform_views_dict.keys()

        self.wf_view_var = tk.StringVar(self._root)
        self.wf_view_var.set(self.waveform_view)  # initial value

        self.om_waveform_view = tk.OptionMenu(self._root, self.wf_view_var,
                                              *tuple(choices),
                                              command=self._wf_views_callback)
        self.om_waveform_view.pack(side=tk.RIGHT)

        l_waveform_view = tk.Label(self._root, text="Waveform")
        l_waveform_view.pack(side=tk.RIGHT)

    def _c_views_callback(self, val):
        self.c_view_var.set(val)
        self.camera_view = val

    def add_tk_camera_views_dropdown(self):
        choices = self.camera_views_dict.keys()

        self.c_view_var = tk.StringVar(self._root)
        self.c_view_var.set(self.camera_view)  # initial value

        self.om_camera_view = tk.OptionMenu(self._root, self.c_view_var,
                                            *tuple(choices),
                                            command=self._c_views_callback)
        self.om_camera_view.pack(side=tk.RIGHT)

        l_camera_view = tk.Label(self._root, text="Camera")
        l_camera_view.pack(side=tk.RIGHT)

    def add_tk_calibration_settings_button(self):
        b = tk.Button(self._root, text="Calibration Parameters",
                      command=self.calibration_popout)
        b.pack(side=tk.RIGHT)

    def calibration_popout(self):
        if self.calib_popout:
            self.calib_popout.destroy()
        self.calib_popout = tk.Toplevel()
        self.calib_popout.title("About this application...")

        if not self.calib_params:
            return

        entry_dict = {}
        i = 0
        d, inv = integrator_dict()
        for key, val in self.calib_params.items():
            tk.Label(self.calib_popout, text=key).grid(row=i)
            entry_dict[key] = tk.Entry(self.calib_popout)
            entry_dict[key].grid(row=i, column=1)
            if val:
                entry_dict[key].delete(0, "end")
                if key == 'integrator':
                    val = inv[val]
                entry_dict[key].insert(0, val)
            i += 1

        def callback():
            cmdline = []
            for param, entry in entry_dict.items():
                value = entry.get()
                if value:
                    cmdkey = '--'+param.replace('integration_', 'integration-')
                    cmdline.append(cmdkey)
                    cmdline.extend(value.split())
            self.set_calibration_params(cmdline)
            self.event = self.event
            self.calib_popout.destroy()

        b = tk.Button(self.calib_popout, text="Calibrate",
                      command=callback)
        b.grid(row=i, column=0, columnspan=2)

    def add_tk_calibrate_checkbox(self):
        var = tk.BooleanVar(self._root)

        def callback():
            self.calibrate = var.get()

        c = tk.Checkbutton(self._root, text="Calibrate", variable=var,
                           command=callback)
        c.pack(side=tk.RIGHT)

    def update_tk_elements(self):
        self.s_event_index.delete(0, "end")
        self.s_event_index.insert(0, self.event_index)
        self.s_event_id.delete(0, "end")
        self.s_event_id.insert(0, self.event_id)
        self.s_telid.delete(0, "end")
        self.s_telid.insert(0, self.telid)

    def update_tk_om(self, om, function, new_choices):
        om['menu'].delete(0, 'end')
        # var = tk.StringVar(self._root)
        # var.set('')
        for choice in new_choices:
            om['menu'].add_command(label=choice,
                                   command=lambda name=choice: function(name))


def main():
    fp = "/Users/Jason/Software/outputs/libCHEC/sky/cameradata_run1594.fits"
    reader = 'targetio'
    # fp = "/Users/Jason/Software/ctapipe/ctapipe-extra/datasets/gamma_test.simtel.gz"
    # reader = 'hessio'
    # fp = "/Users/Jason/Software/outputs/sim_telarray/test_telescope_multitel5.gz"
    # fp = "/Users/Jason/Software/outputs/sim_telarray/meudon_gamma/simtel_runmeudon_gamma_30tel_30deg_1.gz"
    # reader = 'hessio'
    fp = "/Users/Jason/Software/outputs/lab/2016-11-23/run0.fits"
    reader = 'targetio'

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
