#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy import units as u
import argparse
# Ignore warnings about imprecise coordinates (we don't care)
import warnings
warnings.filterwarnings('ignore')

def main():
    parser = argparse.ArgumentParser(description='Create an observing command'\
      ' script for PALFA timing observations.')
    parser.add_argument('--date', '-d', dest='date', metavar='YYYY-MM-DD',\
      required=True, type=str,\
      help='The date (YYYY-MM-DD) of the observation start in AST')
    parser.add_argument('--time', '-t', dest='time', metavar='HH:MM:SS',\
      required=True, type=str,\
      help='The time (HH:MM or HH:MM:SS) of the observation start in AST')
    parser.add_argument('--nhours', '-l', dest='obs_len_hr',\
      required=True, type=float,\
      help='The duration of the observation in decimal hours')
    parser.add_argument('--srcfile', '-s', dest='src_file',\
      required=True, type=str,\
      help='The src file for this observation')
    parser.add_argument('--catfile', '-c', dest='cat_file',\
      required=True, type=str,\
      help='The cat file for this observation')
    parser.add_argument('--outfile', '-o', dest='out_file',\
      required=False, type=str,\
      help='Filename of cmd file to output (cmd file contents printed'\
      ' to screen if this is not used')
    parser.add_argument('--plot', '-p', dest='make_plot',\
      required=False, action='store_true',\
      help='Plot observing path to screen before exiting')
    parser.add_argument('--quiet', '-q', dest='verbose',\
      required=False, action='store_false',\
      help='Suppress most printed output (verbose by default)')
    
    parser.set_defaults(make_plot=False, verbose=True, out_file=None)

    args = parser.parse_args()

    input_date = args.date
    input_time = args.time
    input_obs_len_hr = args.obs_len_hr
    input_src_file = args.src_file
    input_cat_file = args.cat_file

    session = aoSession(input_src_file, input_cat_file, input_date, input_time,
                        input_obs_len_hr)

    session.create_plan(verbose=args.verbose, make_plot=args.make_plot)

    session.make_cmd_file(args.out_file)

# Some things that are useful
za_min = 1.1 # deg
za_max = 19.7 # deg
slew_speed_az = 25. # deg/min
slew_speed_za = 2.5 # deg/min
ao = EarthLocation(x=2390490.0, y=-5564764.0, z=1994727.0, unit='meter')
cal_times = {"search":5, "fold":60}

class aoSource:
    def __init__(self, name, ra_str, dec_str, nsec, mode, receiver,
                 altaz_frame, parfile_path=None, stepsize_sec=None):
        self.name = name
        self.pos = SkyCoord(ra=ra_str, dec=dec_str,
                            unit=('hourangle', 'degree'))
        self.nsec = nsec
        self.mode = mode.lower()
        self.receiver = receiver.lower()
        self.set_altaz_frame(altaz_frame)
        if parfile_path is None:
            self.parfile_path = "/home/gpu/tzpar/%s.par" % name
        else:
            self.parfile_path = parfile_path
        if stepsize_sec is None:
            self.stepsize_sec = (self.altaz.obstime[1] - \
              self.altaz.obstime[0]).jd * 86400.
        else:
            self.stepsize_sec = stepsize_sec

    
    def __repr__(self):
        return "<%s: %s, %s, %s>" % (self.__class__, self.name, self.mode,
                                     self.receiver)

    def set_altaz_frame(self, altaz_frame):
        self.altaz = self.pos.transform_to(altaz_frame)
        
        za = self.get_za()
        times = self.get_times()
        
        self.visible = (za < za_max) * (za > za_min)
        
        changes = np.where(np.diff(self.visible))[0]
        self.ticks_to_up = np.zeros_like(times.jd, dtype=int)
        self.ticks_to_down = np.zeros_like(times.jd, dtype=int)
        if len(changes):
            if not self.visible[0]:
                self.ticks_to_up[0:changes[0]] = np.arange(changes[0], 0, -1)
            else:
                self.ticks_to_down[0:changes[0]] = np.arange(changes[0], 0, -1)
            for ii in range(len(changes)):
                if not self.visible[changes[ii]+1]:
                    if ii < len(changes) - 1:
                        self.ticks_to_up[changes[ii]:changes[ii+1]] = \
                          np.arange(changes[ii+1]-changes[ii], 0, -1)
                    else:
                        self.ticks_to_up[changes[ii]:] = -1
                else:
                    if ii < len(changes) - 1:
                        self.ticks_to_down[changes[ii]:changes[ii+1]] = \
                          np.arange(changes[ii+1]-changes[ii], 0, -1)
                    else:
                        self.ticks_to_down[changes[ii]:] = -1
        else:
            if self.visible.all():
                self.ticks_to_down += -1
            else:
                self.ticks_to_up += -1
        
    def get_za(self):
        return 90. - self.altaz.alt.degree
    
    def get_times(self):
        return self.altaz.obstime
    
    def get_tick_from_time(self, time, absolute_time=False):
        if absolute_time:
            return np.searchsorted(self.get_times().jd, time.jd)
        else:
            tick = int(np.round(time / self.stepsize_sec))
            if tick >= len(self.get_za()):
                return len(self.get_za())-1
            else:
                return tick

    def is_visible(self, time, absolute_time=False):
        """
        If absolute_time, time should be an astropy.time.Time object.
        Otherwise, time should be seconds from the start.
        
        Returns True if visible at this time, False otherwise.
        """
        tick = self.get_tick_from_time(time, absolute_time)
        return self.visible[tick]
        
    def get_seconds_visible(self, time, absolute_time=False):
        """
        If absolute_time, time should be an astropy.time.Time object.
        Otherwise, time should be seconds from the start.
        
        Returns number of seconds before source sets or goes into zenith 
        constraint.
        """
        tick = self.get_tick_from_time(time, absolute_time)
        ticks_to_down = self.ticks_to_down[tick]
        if ticks_to_down >= 0:
            return ticks_to_down*self.stepsize_sec
        else:
            return (len(self.ticks_to_down) - tick)*self.stepsize_sec
    
    def get_total_seconds_remaining(self, time, absolute_time=False):
        """
        If absolute_time, time should be an astropy.time.Time object.
        Otherwise, time should be seconds from the start.        
        
        Unlike get_seconds_visible, this isn't just the time until set *or*
        zenith constraint--it's the total number of seconds left in the whole
        session that this source can be observed.
        """
        tick = self.get_tick_from_time(time, absolute_time)
        nticks_remaining = np.sum(self.visible[tick+1:].astype(int))
        return nticks_remaining * self.stepsize_sec
    
    def get_seconds_not_visible(self, time, absolute_time=False):
        """
        If absolute_time, time should be an astropy.time.Time object.
        Otherwise, time should be seconds from the start.
        
        Returns number of seconds before source rises or comes out of zenith
        constraint.
        """
        tick = self.get_tick_from_time(time, absolute_time)
        ticks_to_up = self.ticks_to_up[tick]
        if ticks_to_up >= 0:
            return ticks_to_up*self.stepsize_sec
        else:
            return (len(self.ticks_to_up) - tick)*self.stepsize_sec
        
    def get_seconds_to_last_seen(self, time, absolute_time=False):
        """
        If absolute_time, time should be an astropy.time.Time object.
        Otherwise, time should be seconds from the start.
        
        Returns number of seconds before source last disappears
        """
        tick = self.get_tick_from_time(time, absolute_time)
        final_tick = len(self.visible) - \
          np.nonzero(np.cumsum(self.visible[::-1]))[0][0] - 1
        return (final_tick - tick) * self.stepsize_sec
            
    def get_seconds_needed(self):
        return self.nsec + cal_times[self.mode]
    
    def slew_time_from(self, other_source, time, absolute_time=False):
        """
        If absolute_time, time should be an astropy.time.Time object.
        Otherwise, time should be seconds from the start.
        
        Returns number of seconds needed to slew to this source from another
        beginning at 'time'.
        """
        tick_now = self.get_tick_from_time(time, absolute_time)
        
        sec_not_visible = self.get_seconds_not_visible(time, absolute_time)
        
        if absolute_time:
            tick = self.get_tick_from_time(time + sec_not_visible*u.second,
                                           absolute_time)
        else:
            tick = self.get_tick_from_time(time + sec_not_visible,
                                           absolute_time)

        other_az = other_source.altaz.az.degree[tick_now]
        other_za = other_source.get_za()[tick_now]

        this_az = self.altaz.az.degree[tick]
        this_za = self.get_za()[tick]

        az_diff = np.abs(other_az - this_az) % 360.
        za_diff = np.abs(other_za - this_za)

        #slew_time_sec = (az_diff/slew_speed_az + za_diff/slew_speed_za) * 60.
        slew_time_sec = max(az_diff/slew_speed_az, za_diff/slew_speed_za) * 60.

        for ii in range(2):
            if absolute_time:
                new_time = time + slew_time_sec*u.second
            else:
                new_time = time + slew_time_sec

            tick = self.get_tick_from_time(new_time, absolute_time)
            this_az = self.altaz.az.degree[tick]
            this_za = self.get_za()[tick]
            az_diff = np.abs(other_az - this_az) % 360.
            za_diff = np.abs(other_za - this_za)
            slew_time_sec = (az_diff/slew_speed_az + za_diff/slew_speed_za) \
              * 60. + sec_not_visible

        #if this_za > za_max or this_za < za_min:
        #    return -1
        #else:
        return slew_time_sec

        
class aoSession:
    """
    srcs_file: name of ASCII file with the following columns:
      (1) source names as they appear in the catalogue file (cat_file)
      (2) number of seconds to observe this source
      (3) "fold" or "search" to choose observing mode
      (4) "LWide" or "430MHz" to choose a receiever
      (5) path to parfile (only needed if mode is "fold")
          (assumed to be /home/gpu/tzpar/[NAME].par if needed and missing)
    cat_file: name of the CIMA catalogue file sources are being read from
      for this session.
    start_date: date of start time (AST) in "YYYY-MM-DD" format 
    start_time: start time (AST) in "HH:MM" or "HH:MM:SS" format
    obslen_hr: total length of session in decimal hours
    stepsize_sec: the altitude and azimuth are calculated every stepsize_sec
      seconds
    """
    def __init__(self, srcs_file, cat_file, start_date, start_time, obslen_hr,
                 stepsize_sec=10):
        self.srcs_file = srcs_file
        self.cat_file = cat_file
        self.stepsize_sec = stepsize_sec
        self.time_steps = np.arange(0, obslen_hr, stepsize_sec/3600.)*u.hour
        self.start_time = Time("%sT%s" % (start_date, start_time), \
          format='isot') + 4.*u.hour
        self.start_time.delta_ut1_utc = 0
        self.obslen_hr = obslen_hr
        self.obslen_sec = int(obslen_hr*3600)
        self.altaz_frame = AltAz(obstime=self.start_time+self.time_steps,\
          location=ao)
        self.plan = []
        
        src_names = []
        src_nsecs = []
        src_modes = []
        src_recvs = []
        src_parfs = []
        with open(srcs_file, 'r') as f:
            for line in f.readlines():
                if line[0] != "#":
                    split_line = line.split()
                    if len(split_line) >= 4:
                        src_name = split_line[0]
                        src_names.append(src_name)
                        src_nsecs.append(int(split_line[1]))
                        mode = split_line[2].lower()
                        if mode in ["search", "fold"]:
                            src_modes.append(mode)
                        else:
                            raise ValueError("Mode input was '%s', must be"\
                              " 'search' or 'fold'." % mode)
                        recv = split_line[3].lower()
                        if recv in ["lwide", "430mhz"]:
                            src_recvs.append(recv)
                        else:
                            raise ValueError("Receiver input was '%s', must"\
                              " be 'LWide' or '430MHz'." % recv)
                        if len(split_line) >= 5:
                            src_parfs.append(split_line[4])
                        else:
                            src_parfs.append("/home/gpu/tzpar/%s.par"%src_name)

        src_ras = {}
        src_decs = {}
        with open(cat_file, 'r') as f:
            for line in f.readlines():
                split_line = line.split()
                if len(split_line) >= 3:
                    cat_name = split_line[0]
                    if cat_name in src_names:
                        src_ras[cat_name] = split_line[1]
                        src_decs[cat_name] = split_line[2]
        
        self.sources = []
        for ii,src in enumerate(src_names):
            self.sources.append(aoSource(src, src_ras[src], src_decs[src], \
              src_nsecs[ii], src_modes[ii], src_recvs[ii], self.altaz_frame, \
              src_parfs[ii], stepsize_sec))
        
    def get_seconds_visible(self, time, absolute_time=False):
        seconds = []
        for ii,src in enumerate(self.sources):
            seconds.append(src.get_seconds_visible(time, absolute_time))
        return np.array(seconds)
    
    def plot_sources(self, time, absolute_time=False, ax=None):
        bgcolor = '0.9'
        gridlines_color = 'white'
        
        if ax is None:
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111)
            ax.set_xlim(25, -25)
            ax.set_ylim(-25, 25)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_axis_bgcolor(bgcolor)
            if absolute_time:
                ax.set_title(time)
            else:
                ax.set_title("%.1f minutes into session" % (time/60.))

        # Let's make a few equatorial coordinate lines
        dec_lines = [-10, 0, 10, 20, 30, 40, 50]
        ra_offs_lines = [-2, -1, 0, 1, 2]
        if not absolute_time:
            actual_time = self.start_time + time*u.second
            actual_time.delta_ut1_utc = 0
        else:
            actual_time = time
        sidereal_time = actual_time.sidereal_time('mean', ao.longitude)
        dec_lines_altaz_frame = AltAz(obstime=actual_time+\
          np.linspace(-3, 3, 100)*u.hour, location=ao)
        ra_lines_altaz_frame = AltAz(obstime=actual_time, location=ao)
        for dec in dec_lines:
            line_coords = SkyCoord(ra=sidereal_time, dec=dec, \
              unit=('hourangle', 'degree'))
            line_altaz = line_coords.transform_to(dec_lines_altaz_frame)
            line_az = line_altaz.az.rad
            line_za = 90. - line_altaz.alt.degree
            line_x = -line_za*np.sin(line_az)
            line_y = line_za*np.cos(line_az)
            ax.plot(line_x, line_y, c=gridlines_color, lw=2, zorder=-10)
        for ra_offs in ra_offs_lines:
            offset_time = actual_time+ra_offs*u.hour
            offset_time.delta_ut1_utc = 0
            sidereal_time = offset_time.sidereal_time('mean', ao.longitude)
            line_coords = SkyCoord(ra=sidereal_time, \
              dec=np.linspace(dec_lines[0], dec_lines[-1], 100), \
              unit=('hourangle', 'degree'))
            line_altaz = line_coords.transform_to(ra_lines_altaz_frame)
            line_az = line_altaz.az.rad
            line_za = 90. - line_altaz.alt.degree
            line_x = -line_za*np.sin(line_az)
            line_y = line_za*np.cos(line_az)
            ax.plot(line_x, line_y, c=gridlines_color, lw=2, zorder=-10)
                
        ax.add_patch(Wedge((0,0), za_max, 0, 360, width=za_max-za_min, \
          facecolor=(1, 1, 1, 0.5)))
            
        src_names = set()
        for src in self.sources:
            if src.name not in src_names:
                tick = src.get_tick_from_time(time, absolute_time)
                az_rad = src.altaz.az.rad
                za = src.get_za()
                xc = -za*np.sin(az_rad)
                yc = za*np.cos(az_rad)
                ax.plot(xc, yc, ':', c='black')
                ax.plot(xc[tick], yc[tick], 'o', c='purple')
                src_names.add(src.name)
                
    def create_plan(self, verbose=True, make_plot=True):
        self.plan = []
        
        done = set()
        t = 0
        az_starts = []
        za_starts = []
        az_ends = []
        za_ends = []
        ra_order = []
        dec_order = []

        ra_hr = np.array([src.pos.ra.hour for src in self.sources])
        st_hr = self.start_time.sidereal_time('mean', ao.longitude).hour
        ha = (st_hr-ra_hr)
        ha[ha > 12] -= 24.
        hour_order = np.argsort(ha)[::-1]
        
        sec_needed = np.array([src.get_seconds_needed() for src in self.sources])
        ra_hours = np.array([src.pos.ra.hour for src in self.sources])
        # assume initial slew is two minutes
        slew_times = np.array([120]*len(self.sources))
        on_source_sec_remaining = sec_needed.sum()
        
        while len(done) < len(self.sources):
            num_left = len(self.sources) - len(done)
            enough_time = True
            
            #if verbose:
            #    print "[%.2f min] Off-source seconds available: %d" % \
            #      (t/60., (self.obslen_sec - t) - on_source_sec_remaining)
            sec_visible = []
            total_sec_remaining = []
            sec_to_last_seen = []
            for ii,src in enumerate(self.sources):
                sec_visible.append(src.get_seconds_visible(t+slew_times[ii]))
                total_sec_remaining.append(src.get_total_seconds_remaining(t+\
                  slew_times[ii]))
                sec_to_last_seen.append(src.get_seconds_to_last_seen(t+\
                  slew_times[ii]))
            sec_visible = np.array(sec_visible)
            total_sec_remaining = np.array(total_sec_remaining)
            sec_to_last_seen = np.array(sec_to_last_seen)
            
            if not total_sec_remaining.any():
                break
                    
            options = []
            for ii in np.where(sec_visible > sec_needed+slew_times)[0]:
                if ii not in done and slew_times[ii] >= 0:
                    options.append(ii)
            if not len(options) and num_left == 1:
                ii_left = set(range(len(self.sources))).difference(done).pop()
                options.append(ii_left)
                enough_time = False
            options = np.array(options)

            if len(options):
                wiggle_room = sec_to_last_seen - sec_needed - slew_times
                this_source = options[np.argmin(wiggle_room[options] * \
                  slew_times[options])]
                if verbose:
                    print "[%.2f min] Of %d option%s, picked %s, %d left"\
                      " (%d sec slew)" % (t/60., len(options), \
                      ['','s'][bool(len(options)-1)], \
                      self.sources[this_source].name, \
                      len(self.sources)-len(done)-1, slew_times[this_source])
                
                self.plan.append(this_source)
                
                t += slew_times[this_source]
                tick = self.sources[this_source].get_tick_from_time(t)
                az_starts.append(self.sources[this_source].altaz.az.degree[tick])
                za_starts.append(self.sources[this_source].get_za()[tick])
                
                t += sec_needed[this_source]
                tick = self.sources[this_source].get_tick_from_time(t)
                az_ends.append(self.sources[this_source].altaz.az.degree[tick])
                za_ends.append(self.sources[this_source].get_za()[tick])
                
                ra_order.append(self.sources[this_source].pos.ra.hour)
                dec_order.append(self.sources[this_source].pos.dec.degree)
                
                on_source_sec_remaining -= sec_needed[this_source]
                
                done.add(this_source)
            else:
                if verbose:
                    print "[%.2f min] Nothing to pick, moving ahead by %d"\
                      " seconds..." % (t/60., self.stepsize_sec)
                t += self.stepsize_sec
            try:
                # This is sloppily placed in a try-except since there's no
                # "this_source" if nothing is found on the first try.  However, 
                # it makes sense for "this_source" to remain the same and for
                # slew times to be updated if nothing is found later.
                # NOTE it would maybe also make sense to add observing time to
                # current pulsar instead of doing nothing, though this could be
                # sensitive to details of the real run
                slew_times = np.array(
                    [src.slew_time_from(self.sources[this_source], t) \
                      for src in self.sources]
                )
            except:
                if verbose:
                    print "Unable to update slew times."
        
        self.plan = np.array(self.plan)

        if verbose:
            print "Total time: %.2f h" % (t / 3600.)
            print "Start time:", self.start_time - 4.*u.hour
            print "End time:", self.start_time - 4.*u.hour + t*u.second
            
        num_missed = len(self.sources) - len(done)
        if num_missed:
            print "Missing %d pulsar%s!" % \
              (num_missed, ['','s'][bool(num_missed-1)])
        if not enough_time:
            print "Not enough time to finish last source!"
        
        if make_plot:
            az_starts = np.array(az_starts)
            za_starts = np.array(za_starts)
            az_ends = np.array(az_ends)
            za_ends = np.array(za_ends)
            ra_order = np.array(ra_order)
            dec_order = np.array(dec_order)

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111)
            ax.set_xlim(25, -25)
            ax.set_ylim(-25, 25)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_axis_bgcolor('0.9')
            self.plot_sources(0, absolute_time=False, ax=ax)
            az_starts_rad = az_starts * np.pi/180.
            az_ends_rad = az_ends * np.pi/180.
            x_starts = -za_starts*np.sin(az_starts_rad)
            y_starts = za_starts*np.cos(az_starts_rad)
            x_ends = -za_ends*np.sin(az_ends_rad)
            y_ends = za_ends*np.cos(az_ends_rad)
            for ii in range(len(x_starts)):
                ax.plot([x_starts[ii], x_ends[ii]], [y_starts[ii], y_ends[ii]],
                        '-', c='blue')
                if ii < len(x_starts)-1:
                    ax.plot([x_ends[ii], x_starts[ii+1]],
                            [y_ends[ii], y_starts[ii+1]], '-', c='orange')
            ax.plot(x_starts, y_starts, 'o', c='blue')
            ax.plot(x_ends, y_ends, 's', c='orange')

            t_start = (self.start_time - 4.*u.hour).datetime
            title_A = t_start.strftime("%Y %b %d %H:%M:%S - ")
            t_end = (self.start_time - 4.*u.hour + t*u.second).datetime
            title_B = t_end.strftime("%H:%M:%S ")
            title_C = "(%.2f hours)" % (t / 3600.)
            ax.set_title(title_A + title_B + title_C)

            plt.show()
        
    def make_cmd_file(self, fname=None):
        """
        If fname=None, just print to screen.
        """
        if not len(self.plan):
            print "No observing plan found--generating one."
            self.create_plan(False, False)
        
        search_lwide = "LOAD PUPPI_LWide_CIMA_Inc_Search.conf\nEXEC change_puppi_scales \"set\" \"04000000\"\n"
        fold_lwide = "LOAD PUPPI_LWide_coherent.conf\nEXEC change_puppi_scales \"set\" \"00800000\"\n"
        fold_430mhz = "LOAD puppi_430MHz_coherent.conf\nEXEC change_puppi_scales \"set\" \"00400000\"\n"
        
        search_setup = "SETUP pulsaron secs = %d loops = 1 caltype = nocal calsecs = %d calmode = on winkcal = off winkcaltype = lcorcal adjpwr = never newfile = one\n"
        fold_setup = "SETUP pulsaron secs = %d loops = 1 caltype = winkcal calsecs = %d calmode = on winkcal = off winkcaltype = lcorcal adjpwr = never newfile = one\n"
        
        cmd_str = ""
        cmd_str += "CATALOG %s\n\n" % self.cat_file
        
        src = self.sources[self.plan[0]]
        current_mode = src.mode
        current_recv = src.receiver
        if current_mode == "search":
            if current_recv == "430mhz":
                raise ValueError("430MHz search mode not supported.")
            elif current_recv == "lwide":
                cmd_str += search_lwide + "\n"
            else:
                raise ValueError("Unknown receiver %s" % current_recv)
            setup_str = search_setup
            #cmd_str += search_setup % (src.nsec, cal_times["search"])
        elif current_mode == "fold":
            if current_recv == "430mhz":
                cmd_str += fold_430mhz + "\n"
            elif current_recv == "lwide":
                cmd_str += fold_lwide + "\n"
            else:
                raise ValueError("Unknown receiver %s" % current_recv)
            setup_str = fold_setup
            #cmd_str += fold_setup % (src.nsec, cal_times["fold"])
        else:
            raise ValueError("Unknown mode %s" % current_mode)
        cmd_str += setup_str % (src.nsec, cal_times[current_mode])
        cmd_str += "SEEK %s\n" % src.name
        if current_mode == "fold":
            cmd_str += "EXEC change_puppi_parfile \"%s\"\n" % (src.parfile_path)
        cmd_str += "ADJUSTPOWER\nADJUSTPOWER\n"
        cmd_str += "EXEC wait_puppi_temporary \"180\" \"Verify PUPPI power levels using 'guppi_adc_hist' in the puppimaster terminal\"\n"
        cmd_str += "PULSARON\n\n"
        
        for ii in range(1, len(self.plan)):
            src = self.sources[self.plan[ii]]
            if src.mode != current_mode or src.receiver != current_recv:
                current_mode = src.mode
                current_recv = src.receiver
                if current_mode == "search":
                    if current_recv == "430mhz":
                        raise ValueError("430MHz search mode not supported.")
                    elif current_recv == "lwide":
                        cmd_str += search_lwide
                    else:
                        raise ValueError("Unknown receiver %s" % current_recv)
                    setup_str = search_setup
                elif current_mode == "fold":
                    if current_recv == "430mhz":
                        cmd_str += fold_430mhz
                    elif current_recv == "lwide":
                        cmd_str += fold_lwide
                    else:
                        raise ValueError("Unknown receiver %s" % current_recv)
                    setup_str = fold_setup
                else:
                    raise ValueError("Unknown mode %s" % current_mode)
                cmd_str += "ADJUSTPOWER\n\n"
                
            cmd_str += "SEEK %s\n" % src.name
            if current_mode == "fold":
                cmd_str += "EXEC change_puppi_parfile \"%s\"\n" % \
                  (src.parfile_path)

            cmd_str += setup_str % (src.nsec, cal_times[current_mode])
            cmd_str += "PULSARON\n\n"
            
        cmd_str += "LOG \"LBand Observation Completed.\"\n"
        
        if fname is None:
            print "# Command file begins below:"
            print ""
            print cmd_str
        else:
            with open(fname, 'w') as f:
                f.write(cmd_str)
            print "Wrote %s" % fname


if __name__ == "__main__":
    main()
