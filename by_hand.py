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

############################################################################################################################################################

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
    
    ############################################################################################################################################################

# The following code will run if this file is executed as a script
if __name__ == "__main__":

    print " "
    
    # Define the time and location of the event
    event_time = obspy.UTCDateTime("2013-05-24T14:56:31")
    event_lat = 52.1357
    event_lon = 151.5688
    event_depth = 632.0
    
    # Get the theoretical arrival times of each phase using TauPy
    model = TauPyModel(model='iasp91')
    arrivals = model.get_travel_times(source_depth_in_km = 632.0, 
                                      distance_in_degree = 82.65, 
                                      phase_list = ["P"])
    p_arrival_time = arrivals[0].time # This time is relative to the event time
    print "P arrival time: %d seconds after the event" % p_arrival_time
    
    # Get the waveforms
    inv, st = get_stations_waveforms("XO", 
                                     "LC15,LD16,LE17,LF18,LG19,LH20,LH22,LI21", 
                                     "BHZ", 
                                     event_time, 
                                     p_arrival_time-60.0, 
                                     p_arrival_time+60.0)
    print "Got %d stations and waveforms." % len(inv[0])

    # Extract the waveform data arrays
    LC15_data = st[0].data
    LD16_data = st[1].data
    LE17_data = st[2].data
    LF18_data = st[3].data
    LG19_data = st[4].data
    LH20_data = st[5].data
    LH22_data = st[6].data
    LI21_data = st[7].data
    
    # Make a time vector
    dt = st[0].stats.delta
    time = [i*dt for i in range(len(LC15_data))]
    
    # Make the traces start 5 seconds before the chosen arrival times
    LC15_data_cut = LC15_data[699:]
    LD16_data_cut = LD16_data[736:]
    LE17_data_cut = LE17_data[774:]
    LF18_data_cut = LF18_data[822:]
    LG19_data_cut = LG19_data[873:]
    LH20_data_cut = LH20_data[906:]
    LH22_data_cut = LH22_data[938:]
    LI21_data_cut = LI21_data[974:]
    
    # Make the traces end 60 seconds after the chosen arrival times
    LC15_data_cut2 = LC15_data_cut[:2600]
    LD16_data_cut2 = LD16_data_cut[:2600]
    LE17_data_cut2 = LE17_data_cut[:2600]
    LF18_data_cut2 = LF18_data_cut[:2600]
    LG19_data_cut2 = LG19_data_cut[:2600]
    LH20_data_cut2 = LH20_data_cut[:2600]
    LH22_data_cut2 = LH22_data_cut[:2600]
    LI21_data_cut2 = LI21_data_cut[:2600]
    
    # Make the traces smaller in amplitude
    amp_factor = 0.2
    LC15_data_cut2[:] *= amp_factor
    LD16_data_cut2[:] *= amp_factor
    LE17_data_cut2[:] *= amp_factor
    LF18_data_cut2[:] *= amp_factor
    LG19_data_cut2[:] *= amp_factor
    LH20_data_cut2[:] *= amp_factor
    LH22_data_cut2[:] *= amp_factor
    LI21_data_cut2[:] *= amp_factor
    
    # Add in each trace's distance
    LC15_data_cut3 = LC15_data_cut2[:] + 75.191335
    LD16_data_cut3 = LD16_data_cut2[:] + 75.349001
    LE17_data_cut3 = LE17_data_cut2[:] + 75.552284
    LF18_data_cut3 = LF18_data_cut2[:] + 75.724926
    LG19_data_cut3 = LG19_data_cut2[:] + 75.935758
    LH20_data_cut3 = LH20_data_cut2[:] + 76.092800
    LH22_data_cut3 = LH22_data_cut2[:] + 76.234206
    LI21_data_cut3 = LI21_data_cut2[:] + 76.374770
    
    # Make a time vector
    dt = st[0].stats.delta
    time = [i*dt for i in range(len(LI21_data_cut2))]
    
    # Plot test
    plt.plot(time, LC15_data_cut3, '-k', time, LD16_data_cut3, '-k', time, LE17_data_cut3, '-k', time, LF18_data_cut3, '-k', time, LG19_data_cut3, '-k', time, LH20_data_cut3, '-k', time, LH22_data_cut3, '-k', time, LI21_data_cut3, '-k')
    plt.xlim(0,65)
    plt.axvline(5)
    plt.xlabel('Time (sec)')
    plt.ylabel('Distance (deg)')
    plt.show()
    
    # Find the intersection times for each slope in the vespagram (and for each distance at each slope)
    slopes = [-1.0, -4.0/3.0, -2.0, -4.0, 4.0, 2.0, 4.0/3.0, 1.0]
    distances = [75.191335, 75.349001, 75.552284, 75.724926, 75.935758, 76.092800, 76.234206, 76.374770]
    for slope in slopes:
        for distance in distances:
            intersection_time = 5.0 - (75.724926 - distance)/(slope)
            intersection_index = int(round(intersection_time * (1.0/dt)))
            #print "For the distance of %f and the slope of %f, the intersection index is %f " % (distance, slope, intersection_index)
     
    # Make the summed waveform for the first slope: -1.0
    LC15_data_slope_01 = LC15_data_cut2[47:]
    LD16_data_slope_01 = LD16_data_cut2[41:]
    LE17_data_slope_01 = LE17_data_cut2[33:]
    LF18_data_slope_01 = LF18_data_cut2[26:]
    LG19_data_slope_01 = LG19_data_cut2[18:]
    LH20_data_slope_01 = LH20_data_cut2[11:]
    LH22_data_slope_01 = LH22_data_cut2[6:]
    LI21_data_slope_01 = LI21_data_cut2[0:]
    LC15_data_slope_01 = LC15_data_slope_01[:len(LC15_data_slope_01)]
    LD16_data_slope_01 = LD16_data_slope_01[:len(LC15_data_slope_01)]
    LE17_data_slope_01 = LE17_data_slope_01[:len(LC15_data_slope_01)]
    LF18_data_slope_01 = LF18_data_slope_01[:len(LC15_data_slope_01)]
    LG19_data_slope_01 = LG19_data_slope_01[:len(LC15_data_slope_01)]
    LH20_data_slope_01 = LH20_data_slope_01[:len(LC15_data_slope_01)]
    LH22_data_slope_01 = LH22_data_slope_01[:len(LC15_data_slope_01)]
    LI21_data_slope_01 = LI21_data_slope_01[:len(LC15_data_slope_01)]
    summed_trace_for_slope_01 = [sum(x) for x in zip(LC15_data_slope_01, LD16_data_slope_01, LE17_data_slope_01, LF18_data_slope_01, LG19_data_slope_01, LH20_data_slope_01, LH22_data_slope_01, LI21_data_slope_01)]
    summed_trace_for_slope_01 = np.array(summed_trace_for_slope_01)
    #plt.plot(summed_trace_for_slope_01)
    #plt.show()
    
    # Make the summed waveform for the second slope: -1.333333
    LC15_data_slope_02 = LC15_data_cut2[35:]
    LD16_data_slope_02 = LD16_data_cut2[30:]
    LE17_data_slope_02 = LE17_data_cut2[24:]
    LF18_data_slope_02 = LF18_data_cut2[19:]
    LG19_data_slope_02 = LG19_data_cut2[13:]
    LH20_data_slope_02 = LH20_data_cut2[8:]
    LH22_data_slope_02 = LH22_data_cut2[4:]
    LI21_data_slope_02 = LI21_data_cut2[0:]
    LC15_data_slope_02 = LC15_data_slope_02[:len(LC15_data_slope_02)]
    LD16_data_slope_02 = LD16_data_slope_02[:len(LC15_data_slope_02)]
    LE17_data_slope_02 = LE17_data_slope_02[:len(LC15_data_slope_02)]
    LF18_data_slope_02 = LF18_data_slope_02[:len(LC15_data_slope_02)]
    LG19_data_slope_02 = LG19_data_slope_02[:len(LC15_data_slope_02)]
    LH20_data_slope_02 = LH20_data_slope_02[:len(LC15_data_slope_02)]
    LH22_data_slope_02 = LH22_data_slope_02[:len(LC15_data_slope_02)]
    LI21_data_slope_02 = LI21_data_slope_02[:len(LC15_data_slope_02)]
    summed_trace_for_slope_02 = [sum(x) for x in zip(LC15_data_slope_02, LD16_data_slope_02, LE17_data_slope_02, LF18_data_slope_02, LG19_data_slope_02, LH20_data_slope_02, LH22_data_slope_02, LI21_data_slope_02)]
    summed_trace_for_slope_02 = np.array(summed_trace_for_slope_02)
    #plt.plot(summed_trace_for_slope_02)
    #plt.show()
    
    # Make the summed waveform for the third slope: -2.0
    LC15_data_slope_03 = LC15_data_cut2[24:]
    LD16_data_slope_03 = LD16_data_cut2[21:]
    LE17_data_slope_03 = LE17_data_cut2[16:]
    LF18_data_slope_03 = LF18_data_cut2[13:]
    LG19_data_slope_03 = LG19_data_cut2[9:]
    LH20_data_slope_03 = LH20_data_cut2[6:]
    LH22_data_slope_03 = LH22_data_cut2[3:]
    LI21_data_slope_03 = LI21_data_cut2[0:]
    LC15_data_slope_03 = LC15_data_slope_03[:len(LC15_data_slope_03)]
    LD16_data_slope_03 = LD16_data_slope_03[:len(LC15_data_slope_03)]
    LE17_data_slope_03 = LE17_data_slope_03[:len(LC15_data_slope_03)]
    LF18_data_slope_03 = LF18_data_slope_03[:len(LC15_data_slope_03)]
    LG19_data_slope_03 = LG19_data_slope_03[:len(LC15_data_slope_03)]
    LH20_data_slope_03 = LH20_data_slope_03[:len(LC15_data_slope_03)]
    LH22_data_slope_03 = LH22_data_slope_03[:len(LC15_data_slope_03)]
    LI21_data_slope_03 = LI21_data_slope_03[:len(LC15_data_slope_03)]
    summed_trace_for_slope_03 = [sum(x) for x in zip(LC15_data_slope_03, LD16_data_slope_03, LE17_data_slope_03, LF18_data_slope_03, LG19_data_slope_03, LH20_data_slope_03, LH22_data_slope_03, LI21_data_slope_03)]
    summed_trace_for_slope_03 = np.array(summed_trace_for_slope_03)
    #plt.plot(summed_trace_for_slope_03)
    #plt.show()
    
    # Make the summed waveform for the fourth slope: -4.0
    LC15_data_slope_04 = LC15_data_cut2[11:]
    LD16_data_slope_04 = LD16_data_cut2[10:]
    LE17_data_slope_04 = LE17_data_cut2[8:]
    LF18_data_slope_04 = LF18_data_cut2[6:]
    LG19_data_slope_04 = LG19_data_cut2[4:]
    LH20_data_slope_04 = LH20_data_cut2[2:]
    LH22_data_slope_04 = LH22_data_cut2[1:]
    LI21_data_slope_04 = LI21_data_cut2[0:]
    LC15_data_slope_04 = LC15_data_slope_04[:len(LC15_data_slope_04)]
    LD16_data_slope_04 = LD16_data_slope_04[:len(LC15_data_slope_04)]
    LE17_data_slope_04 = LE17_data_slope_04[:len(LC15_data_slope_04)]
    LF18_data_slope_04 = LF18_data_slope_04[:len(LC15_data_slope_04)]
    LG19_data_slope_04 = LG19_data_slope_04[:len(LC15_data_slope_04)]
    LH20_data_slope_04 = LH20_data_slope_04[:len(LC15_data_slope_04)]
    LH22_data_slope_04 = LH22_data_slope_04[:len(LC15_data_slope_04)]
    LI21_data_slope_04 = LI21_data_slope_04[:len(LC15_data_slope_04)]
    summed_trace_for_slope_04 = [sum(x) for x in zip(LC15_data_slope_04, LD16_data_slope_04, LE17_data_slope_04, LF18_data_slope_04, LG19_data_slope_04, LH20_data_slope_04, LH22_data_slope_04, LI21_data_slope_04)]
    summed_trace_for_slope_04 = np.array(summed_trace_for_slope_04)
    #plt.plot(summed_trace_for_slope_04)
    #plt.show()
    
    # Make the summed waveform for the center case: the initial plot w/ the vertical line
    summed_trace_for_center = [sum(x) for x in zip(LC15_data_cut2, LD16_data_cut2, LE17_data_cut2, LF18_data_cut2, LG19_data_cut2, LH20_data_cut2, LH22_data_cut2, LI21_data_cut2)]
    summed_trace_for_center = np.array(summed_trace_for_center)
    #plt.plot(summed_trace_for_center)
    #plt.show()
    
    # Make the summed waveform for the fifth slope: 4.0
    LC15_data_slope_05 = LC15_data_cut2[0:]
    LD16_data_slope_05 = LD16_data_cut2[1:]
    LE17_data_slope_05 = LE17_data_cut2[3:]
    LF18_data_slope_05 = LF18_data_cut2[5:]
    LG19_data_slope_05 = LG19_data_cut2[7:]
    LH20_data_slope_05 = LH20_data_cut2[9:]
    LH22_data_slope_05 = LH22_data_cut2[10:]
    LI21_data_slope_05 = LI21_data_cut2[11:]
    LC15_data_slope_05 = LC15_data_slope_05[:len(LI21_data_slope_05)]
    LD16_data_slope_05 = LD16_data_slope_05[:len(LI21_data_slope_05)]
    LE17_data_slope_05 = LE17_data_slope_05[:len(LI21_data_slope_05)]
    LF18_data_slope_05 = LF18_data_slope_05[:len(LI21_data_slope_05)]
    LG19_data_slope_05 = LG19_data_slope_05[:len(LI21_data_slope_05)]
    LH20_data_slope_05 = LH20_data_slope_05[:len(LI21_data_slope_05)]
    LH22_data_slope_05 = LH22_data_slope_05[:len(LI21_data_slope_05)]
    LI21_data_slope_05 = LI21_data_slope_05[:len(LI21_data_slope_05)]
    summed_trace_for_slope_05 = [sum(x) for x in zip(LC15_data_slope_05, LD16_data_slope_05, LE17_data_slope_05, LF18_data_slope_05, LG19_data_slope_05, LH20_data_slope_05, LH22_data_slope_05, LI21_data_slope_05)]
    summed_trace_for_slope_05 = np.array(summed_trace_for_slope_05)
    #plt.plot(summed_trace_for_slope_05)
    #plt.show()
    
    # Make the summed waveform for the sixth slope: 2.0
    LC15_data_slope_06 = LC15_data_cut2[0:]
    LD16_data_slope_06 = LD16_data_cut2[3:]
    LE17_data_slope_06 = LE17_data_cut2[8:]
    LF18_data_slope_06 = LF18_data_cut2[11:]
    LG19_data_slope_06 = LG19_data_cut2[15:]
    LH20_data_slope_06 = LH20_data_cut2[18:]
    LH22_data_slope_06 = LH22_data_cut2[21:]
    LI21_data_slope_06 = LI21_data_cut2[24:]
    LC15_data_slope_06 = LC15_data_slope_06[:len(LI21_data_slope_06)]
    LD16_data_slope_06 = LD16_data_slope_06[:len(LI21_data_slope_06)]
    LE17_data_slope_06 = LE17_data_slope_06[:len(LI21_data_slope_06)]
    LF18_data_slope_06 = LF18_data_slope_06[:len(LI21_data_slope_06)]
    LG19_data_slope_06 = LG19_data_slope_06[:len(LI21_data_slope_06)]
    LH20_data_slope_06 = LH20_data_slope_06[:len(LI21_data_slope_06)]
    LH22_data_slope_06 = LH22_data_slope_06[:len(LI21_data_slope_06)]
    LI21_data_slope_06 = LI21_data_slope_06[:len(LI21_data_slope_06)]
    summed_trace_for_slope_06 = [sum(x) for x in zip(LC15_data_slope_06, LD16_data_slope_06, LE17_data_slope_06, LF18_data_slope_06, LG19_data_slope_06, LH20_data_slope_06, LH22_data_slope_06, LI21_data_slope_06)]
    summed_trace_for_slope_06 = np.array(summed_trace_for_slope_06)
    #plt.plot(summed_trace_for_slope_06)
    #plt.show()
    
    # Make the summed waveform for the seventh slope: 1.333333
    LC15_data_slope_07 = LC15_data_cut2[0:]
    LD16_data_slope_07 = LD16_data_cut2[5:]
    LE17_data_slope_07 = LE17_data_cut2[11:]
    LF18_data_slope_07 = LF18_data_cut2[16:]
    LG19_data_slope_07 = LG19_data_cut2[22:]
    LH20_data_slope_07 = LH20_data_cut2[27:]
    LH22_data_slope_07 = LH22_data_cut2[31:]
    LI21_data_slope_07 = LI21_data_cut2[35:]
    LC15_data_slope_07 = LC15_data_slope_07[:len(LI21_data_slope_07)]
    LD16_data_slope_07 = LD16_data_slope_07[:len(LI21_data_slope_07)]
    LE17_data_slope_07 = LE17_data_slope_07[:len(LI21_data_slope_07)]
    LF18_data_slope_07 = LF18_data_slope_07[:len(LI21_data_slope_07)]
    LG19_data_slope_07 = LG19_data_slope_07[:len(LI21_data_slope_07)]
    LH20_data_slope_07 = LH20_data_slope_07[:len(LI21_data_slope_07)]
    LH22_data_slope_07 = LH22_data_slope_07[:len(LI21_data_slope_07)]
    LI21_data_slope_07 = LI21_data_slope_07[:len(LI21_data_slope_07)]
    summed_trace_for_slope_07 = [sum(x) for x in zip(LC15_data_slope_07, LD16_data_slope_07, LE17_data_slope_07, LF18_data_slope_07, LG19_data_slope_07, LH20_data_slope_07, LH22_data_slope_07, LI21_data_slope_07)]
    summed_trace_for_slope_07 = np.array(summed_trace_for_slope_07)
    #plt.plot(summed_trace_for_slope_07)
    #plt.show()
    
    # Make the summed waveform for the eighth slope: 1.0
    LC15_data_slope_08 = LC15_data_cut2[0:]
    LD16_data_slope_08 = LD16_data_cut2[6:]
    LE17_data_slope_08 = LE17_data_cut2[14:]
    LF18_data_slope_08 = LF18_data_cut2[21:]
    LG19_data_slope_08 = LG19_data_cut2[29:]
    LH20_data_slope_08 = LH20_data_cut2[36:]
    LH22_data_slope_08 = LH22_data_cut2[41:]
    LI21_data_slope_08 = LI21_data_cut2[47:]
    LC15_data_slope_08 = LC15_data_slope_08[:len(LI21_data_slope_08)]
    LD16_data_slope_08 = LD16_data_slope_08[:len(LI21_data_slope_08)]
    LE17_data_slope_08 = LE17_data_slope_08[:len(LI21_data_slope_08)]
    LF18_data_slope_08 = LF18_data_slope_08[:len(LI21_data_slope_08)]
    LG19_data_slope_08 = LG19_data_slope_08[:len(LI21_data_slope_08)]
    LH20_data_slope_08 = LH20_data_slope_08[:len(LI21_data_slope_08)]
    LH22_data_slope_08 = LH22_data_slope_08[:len(LI21_data_slope_08)]
    LI21_data_slope_08 = LI21_data_slope_08[:len(LI21_data_slope_08)]
    summed_trace_for_slope_08 = [sum(x) for x in zip(LC15_data_slope_08, LD16_data_slope_08, LE17_data_slope_08, LF18_data_slope_08, LG19_data_slope_08, LH20_data_slope_08, LH22_data_slope_08, LI21_data_slope_08)]
    summed_trace_for_slope_08 = np.array(summed_trace_for_slope_08)
    #plt.plot(summed_trace_for_slope_08)
    #plt.show()
    
    #print len(summed_trace_for_slope_01)
    #print len(summed_trace_for_slope_02)
    #print len(summed_trace_for_slope_03)
    #print len(summed_trace_for_slope_04)
    #print len(summed_trace_for_center)
    #print len(summed_trace_for_slope_05)
    #print len(summed_trace_for_slope_06)
    #print len(summed_trace_for_slope_07)
    #print len(summed_trace_for_slope_08)
    
    # Downsize each summed trace's amplitude
    amp_factor_2 = 0.2
    summed_trace_for_slope_01 *= amp_factor_2
    summed_trace_for_slope_02 *= amp_factor_2
    summed_trace_for_slope_03 *= amp_factor_2
    summed_trace_for_slope_04 *= amp_factor_2
    summed_trace_for_center *= amp_factor_2
    summed_trace_for_slope_05 *= amp_factor_2
    summed_trace_for_slope_06 *= amp_factor_2
    summed_trace_for_slope_07 *= amp_factor_2
    summed_trace_for_slope_08 *= amp_factor_2
    
    # Add/subtract the necessary slowness (1/slope) from each summed trace's values
    summed_trace_for_slope_01 += -1.0
    summed_trace_for_slope_02 += -0.75
    summed_trace_for_slope_03 += -0.5
    summed_trace_for_slope_04 += -0.25
    summed_trace_for_center += 0.0
    summed_trace_for_slope_05 += 0.25
    summed_trace_for_slope_06 += 0.50
    summed_trace_for_slope_07 += 0.75
    summed_trace_for_slope_08 += 1.0
    
    # Trim each summed trace to the length of the shortest one
    sum01 = summed_trace_for_slope_01[:2553]
    sum02 = summed_trace_for_slope_02[:2553]
    sum03 = summed_trace_for_slope_03[:2553]
    sum04 = summed_trace_for_slope_04[:2553]
    sum_center = summed_trace_for_center[:2553]
    sum05 = summed_trace_for_slope_05[:2553]
    sum06 = summed_trace_for_slope_06[:2553]
    sum07 = summed_trace_for_slope_07[:2553]
    sum08 = summed_trace_for_slope_08[:2553]
    
    # Make a time vector
    time_for_vespagram = [i*dt for i in range(len(sum01))]
    
    # Plot the vespagram
    plt.plot(time_for_vespagram, sum01, '-k', time_for_vespagram, sum02, '-k', time_for_vespagram, sum03, '-k', time_for_vespagram, sum04, '-k', time_for_vespagram, sum_center, '-k', time_for_vespagram, sum05, '-k', time_for_vespagram, sum06, '-k', time_for_vespagram, sum07, '-k', time_for_vespagram, sum08, '-k')
    plt.xlabel('Time (sec)')
    plt.ylabel('Slowness (sec/deg)')
    plt.show()
    
    # Plot the vespagram as a contour plot
    #slownesses = [-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
    #matrix_of_summed_waveforms = [sum01, sum02, sum03, sum04, sum_center, sum05, sum06, sum07, sum08]
    #trace_matrix = np.zeros((len(slownesses), len(sum01)))
    #X, Y = np.meshgrid(time_for_vespagram, slownesses)
    #for i in range(len(slownesses)):
        #trace_matrix[i,:] = matrix_of_summed_waveforms[i]

    # Create the plot itself
    #plt.figure()
    #matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    #plt.contour(X, Y, trace_matrix, 10, colors='k')
    #plt.xlim(min(time_for_vespagram), max(time_for_vespagram))
    #plt.ylim(min(slownesses), max(slownesses))
    #plt.xlabel('Time (sec)')
    #plt.ylabel('Slowness (sec/deg)')
    #plt.show()
