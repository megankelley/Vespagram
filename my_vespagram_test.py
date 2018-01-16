# Created by M. Kelley on 20 November 2017
# Output: a linear contour vespagram

import obspy
from obspy.clients.fdsn import Client
from obspy.core.inventory import read_inventory
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from obspy.geodetics.base import locations2degrees
import numpy as np
import matplotlib
from obspy.taup import TauPyModel

#######################################

def get_stations_waveforms(network, stations, channel, event_time, waveform_start, waveform_end):
    
    """
    Retrieve an inventory of stations and a stream of waveforms from IRIS
    
    Args:
        network:        String of ONE network
        stations:       Comma-separated string of station codes in network
        channel:        String of channel name
        event_time:     obspy.UTCDateTime of event
        waveform_start: Integer seconds after event_time to start waveform
        waveform_end:   Integer seconds after event_time to end waveform
       
    Returns: obspy.Inventory and obspy.Stream for a chosen event
    """
    
    # Retrieve the information
    c = Client("IRIS")
    inv = c.get_stations(network=network, station=stations, channel=channel, level="response")
    st = c.get_waveforms(network=network, station=stations, location="", channel=channel, starttime=event_time+waveform_start, endtime=event_time+waveform_end) # Units of seconds

    # Filter, remove the instrument response, normalize the data values to a maximum value of 1
    st.attach_response(inv)
    st.detrend("linear")
    pre_filt = [0.005, 0.04, 5.0, 6.0]
    st.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP")
    st.normalize()
    
    return inv, st

#######################################
    
def make_record_section(stream, network, event_lat, event_lon, amp_factor):

    """
    Makes a plot of data traces arranged by their distance from an event
    
    Args:
        stream:         obspy.Stream object, containing the waveform data
        network:        String of ONE network
        event_lat:      Latitude of the event, decimal format
        event_lon:      Longitude of the event, decimal format
        amp_factor:     A double format amplification factor for the amplitude of the traces in the plot
       
    Returns: a matplotlib figure saved as a .png file
    """
    
    # Clear any existing figures
    plt.clf()
    
    # Make the record section one trace at a time
    for station, trace in zip(network, stream):
        dist = locations2degrees(event_lat, event_lon, station.latitude, station.longitude)
        data = trace.data
        data *= amp_factor
        data[:] += dist
        dt = stream[0].stats.delta
        time = [i*dt for i in range(len(data))]
        plt.plot(time, data, '-k', label=station.code)
        
    plt.xlabel('Time (sec)')
    plt.ylabel('Distance (deg)')
    plt.show()
    
    # Save the figure
    plt.savefig('record_section.png')

#######################################

def get_slowness_intercept_pairs(stream, network, event_lat, event_lon, event_depth, delta_upper, delta_lower, num_samples):
    """
    Calculates the intercepts and slownesses of lines rotated around the central trace's P arrival time
    
    Args:
        stream:         obspy.Stream object, containing the waveform data
        network:        String of ONE network
        event_lat:      Latitude of the event, decimal format
        event_lon:      Longitude of the event, decimal format
        event_depth:    Depth of the event, decimal format, units of km
        delta_upper:    Range of slownesses above average
        delta_lower:    Range of slownesses below average
        num_samples:    Number of slowness samples
    
    Returns: a list of (slope, intercept) pairs for each slowness value
    """
    
    # Preliminary information
    distances = []
    p_times = []
    dt = stream[0].stats.delta
    model = TauPyModel(model='iasp91')
    
    # Calculate the distance and P arrival time for each station
    for station, trace in zip(network, stream):
        dist = locations2degrees(event_lat, event_lon, station.latitude, station.longitude)
        distances.append(dist)
        arrival = model.get_travel_times(event_depth, dist, ["P"])
        p_times.append(arrival[0].time) # [0] corresponds to the first arrival in the list, the P wave
        
    # Convert these P arrival times from seconds after the event to seconds after the waveform start
    for time, trace, i in zip(p_times,stream, range(len(p_times))):
        p_times[i] = time - (trace.stats.starttime - event_time)
    
    # Find the nearest station to the event, calculate its P arrival time
    nearest_station_index = np.argmin(distances)
    nearest_dist = distances[nearest_station_index]
    nearest_p_time = p_times[nearest_station_index]
    #print "P time of nearest station: %f seconds" % nearest_p_time
    #print "Distance of nearest station: %f degrees" % nearest_dist
    
    # Find the furthest station from the event, calculate its P arrival time
    furthest_station_index = np.argmax(distances)
    furthest_dist = distances[furthest_station_index]
    furthest_p_time = p_times[furthest_station_index]
    #print "P time of furtest station: %f seconds" % furthest_p_time
    #print "Distance of furthest station: %f degrees" % furthest_dist
    
    # Choose the station closest to the center of the distance range, save its distance and p arrival time
    avg_dist = float(sum(distances))/len(distances)
    ref_distance = min(distances, key=lambda x:abs(x-avg_dist))
    ref_time = p_times[distances.index(ref_distance)]
    
    # Calculate the slope and intercept of average slowness
    velocity = (furthest_dist - nearest_dist)/(furthest_p_time - nearest_p_time)
    avg_p_slowness = 1.0/velocity
    avg_intercept = ref_distance - velocity*ref_time
    #print "Average P slowness, from inside the function: %f seconds/degree" % avg_p_slowness
    #print "Intercept for average p slowness: %f" % avg_intercept
    
    # Define a vector of slowness values to sum over
    slownesses = np.linspace(avg_p_slowness-delta_lower, avg_p_slowness+delta_upper, num_samples)
    
    # For each of the slowness values, calculate intercept
    slowness_intercepts = []
    for slowness in slownesses:
        # b' = t*(m - m') + b where t is the center trace's p arrival time, m's are slopes, and b's are intercepts
        intercept = ref_time*(1.0/avg_p_slowness - 1.0/slowness) + avg_intercept
        slowness_intercepts.append((slowness, intercept))
    
    return slowness_intercepts

#######################################

def stack_waveforms(stream, network, event_lat, event_lon, slowness, intercept):

    """
    Shifts a group of waveforms according to a given slowness and sums them
    
    Args:
        stream:         obspy.Stream object, containing the waveform data
        network:        String of ONE network
        event_lat:      Latitude of the event, decimal format
        event_lon:      Longitude of the event, decimal format
        slowness:       The slowness (1/velocity = 1/slope) of the waveform across the array in sec/degree
        intercept:      The y-intercept of the line that represents the arriving wavefront
       
    Returns: one summed waveform
    """
    
    # Empty lists of things I need to make the vespagram later
    names = []
    data_shifts = []
    
    dt = stream[0].stats.delta
    
    # Fill the lists
    for station, trace in zip(network, stream):
        dist = locations2degrees(event_lat, event_lon, station.latitude, station.longitude)
        cross_pt = (dist-intercept)*slowness
        slice_index = int(cross_pt/dt)
        if len(trace.data)-slice_index < 1000:
            print "%s doesn't have enough samples (length %d, slice_index %d)" % (station.code, len(trace.data), slice_index)
        else:
            names.append(station.code)
            data_shifts.append(trace.data[slice_index:])
        
    # Trim data to the same size and reduce
    short_trace_len = min([len(x) for x in data_shifts])
    traces_sum = np.zeros(short_trace_len)
    #print("short_trace_len: {}".format(short_trace_len))
    
    for data in data_shifts:
        data_trim = data[:short_trace_len]
        traces_sum += data_trim
    
    inf_norm = np.max(traces_sum)
    traces_sum[:] /= inf_norm
    
    return traces_sum

#######################################

# The following code will run if this file is executed as a script
if __name__ == "__main__":

    print " "
    
    # Define the time and location of the event
    event_time = obspy.UTCDateTime("2013-05-24T14:56:31")
    event_lat = 52.1357
    event_lon = 151.5688
    event_depth = 632.0
    
    # Get the theoretical arrival times of each phase using TauPy
    model = TauPyModel(model='prem')
    arrivals = model.get_travel_times(source_depth_in_km = 632.0, distance_in_degree = 82.65, phase_list = ["P", "PP", "PcP"])
    p_arrival_time = arrivals[0].time
    
    # Test get_stations_waveforms function
    inv, st = get_stations_waveforms("Z9", "W08,W09,W10,W11,W13,W14,W15A,W16,W18,W19", "BHZ", event_time, p_arrival_time-30.0, p_arrival_time+60.0)
    print "Got %d stations and waveforms." % len(inv[0])

    # Test the make_record_section function
    #make_record_section(st, inv[0], event_lat, event_lon, 0.1)
    #print "Made record section."
    
    # Define a range of slopes (1/slownesses) to investigate
    #slownesses = np.linspace(5.5, 6.5, 40)
    #slopes = [1.0/slowness for slowness in slownesses]
    
    slowness_intercepts = get_slowness_intercept_pairs(st, inv[0], event_lat, event_lon, event_depth, 2.0, 0.5, 20)
    print slowness_intercepts

    # Use the stack_waveforms function to make a summed trace for each of the defined slopes
    trace_sums = []
    for slowness, intercept in slowness_intercepts:
        trace_sums.append(stack_waveforms(st, inv[0], event_lat, event_lon, slowness, intercept))
    print "Summed the shifted traces at each of the different slopes/slownesses."
    
    # Find the length of the shortest trace
    short_trace_len = min([len(x) for x in trace_sums])
    
    # Trim each of the traces in trace_sums to be the length of the shortest trace
    for i, trace in enumerate(trace_sums):
        trace_sums[i] = trace[:short_trace_len]
    print "Trimmed all the traces to the length of the shortest trace."
    
    # Axes for plots
    dt = st[0].stats.delta
    time = [i*dt for i in range(short_trace_len)]
    slownesses = [si[0] for si in slowness_intercepts]
    
    # Test slowness intercepts plot
    #plt.clf()
    #for slow_int, trace in zip(slowness_intercepts, trace_sums):
    #    slope = 1.0/slow_int[0]
    #    yint = slow_int[1]
    #    y = slope*np.asarray(time) + yint
    #    plt.plot(time, y, '-r')
    #for station, trace in zip(inv[0], st):
    #    dist = locations2degrees(event_lat, event_lon, station.latitude, station.longitude)
    #    data = trace.data
    #    data *= 1.0
    #    data[:] += dist
    #    dt = st[0].stats.delta
    #    t = [i*dt for i in range(len(data))]
    #    plt.plot(t, data, '-k')
    #plt.xlabel('Time (sec)')
    #plt.ylabel('Distance (deg)')
    #plt.show()
    
    # Make a wiggle plot
    plt.clf()
    amp_factor = slownesses[1]-slownesses[0]
    amp_factor = 0.05
    for slowness, trace in zip(slownesses, trace_sums):
        trace = trace * amp_factor 
        trace = np.array(trace)
        for i in range(len(trace)):
            trace[i] += slowness
        plt.plot(time, trace, '-k')
    
    plt.xlim(min(time), max(time))
    plt.ylim(min(slownesses)-0.1, max(slownesses)+0.1)
    plt.xlabel('Time (sec)')
    plt.ylabel('Slowness (sec/deg)')
    plt.show()
    
    # Make a contour plot
    slownesses = [si[0] for si in slowness_intercepts]
    trace_matrix = np.zeros((len(slownesses), short_trace_len))
    X, Y = np.meshgrid(time, slownesses)
    for i in range(len(slownesses)):
        trace_matrix[i,:] = trace_sums[i]
    
    plt.figure()
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    plt.contour(X, Y, trace_matrix, 10, colors='k')
    plt.xlim(min(time), max(time))
    plt.ylim(min(slownesses), max(slownesses))
    plt.xlabel('Time (sec)')
    plt.ylabel('Slowness (sec/deg)')
    plt.show()
    
    
