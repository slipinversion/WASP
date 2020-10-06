# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 17:34:16 2017

@author: pk
"""

# Manual data acquisition


import time
import queue
import threading
from obspy.clients.fdsn import Client
from obspy.clients.iris import Client as IRIS_Client
from obspy.core.inventory.inventory import Inventory
import os


def acquisition(event_time, lat_ep, lon_ep, depth, data_to_use, eew=False):
    """
    """
    t1 = event_time - 3 * 60
    t2 = event_time + 81 * 60
    t3 = event_time - 60 if not eew else event_time - 5 * 60
    t4 = event_time + 5 * 60 if not eew else event_time + 140 * 60
    client_iris = Client("IRIS")
    try:
        client_gfz = Client("GFZ") #For stations belonging to CX network
    except:
        pass
    inventory = Inventory([], None)
    if 'strong' in data_to_use:
        networks = "C,C1,II,IU"
        inventory = client_iris.get_stations(
            starttime=event_time - 60, endtime=event_time + 5*60,
            network=networks, channel="HN*", level="response",
            maxradius=10, latitude=lat_ep, longitude=lon_ep)
        networks = "CX"
        try:
            inventory = inventory + client_gfz.get_stations(
                starttime=event_time - 60, endtime=event_time + 5*60,
                network=networks, channel="HL*", level="response",
                maxradius=10, latitude=lat_ep, longitude=lon_ep)
            print('we have the inventory')
        except:
            pass
    inventory_tele = Inventory([], None)
    if 'tele' in data_to_use:
        networks = "II,G,IU,GE"
        max_dist = 90
        inventory_tele = client_iris.get_stations(
            starttime=event_time - 10*60, endtime=event_time + 120*60,
            network=networks, channel="BH*", level="response", minradius=30,
            maxradius=max_dist, latitude=lat_ep, longitude=lon_ep)
    
    cola_tr = queue.Queue()
    cola_mt = queue.Queue()
    
    iris_client = IRIS_Client()
    for network in inventory_tele:
        netwk = network.code
        for station in network:
            statn = station.code
            for canal in station:
                loc_code = canal.location_code
                channel = canal.code
                sac_dict = __get_channel_information_manual(
                    canal, lat_ep, lon_ep, depth)
                t_tr = threading.Thread(
                    target=lambda q, a, b, c, d, e, f, g, h: q.put(
                        worker2(a, b, c, d, e, f, g, h)),
                    args = (cola_tr, client_iris, netwk, statn, loc_code,
                            channel, sac_dict, t1, t2))
                t_mt = threading.Thread(
                    target=lambda q, a, b, c, d, e, f, g: q.put(
                        worker3(a, b, c, d, e, f, g)),
                    args = (cola_mt, iris_client, netwk, statn, loc_code,
                            channel, t1, t2))
                t_tr.daemon = True
                t_tr.start()
                t_mt.daemon = True
                t_mt.start()
    time.sleep(100)

    cola_tr = queue.Queue()
    cola_mt = queue.Queue()

    for network in inventory:
        netwk = network.code
        for station in network:
            statn = station.code
            if statn in ['AC02']:
                continue
            for canal in station:
                loc_code = canal.location_code
                channel = canal.code
                sac_dict = __get_channel_information_manual(
                    canal, lat_ep, lon_ep, depth)
#                worker1(netwk, statn, loc_code, channel, sac_dict, t3, t4)
#                worker4(canal, netwk, statn, loc_code, channel, t3, t4)
                if not netwk=='CX':
                    t_tr = threading.Thread(
                        target=lambda q, a, b, c, d, e, f, g, h:\
                            q.put(worker2(a, b, c, d, e, f, g, h)),
                        args = (cola_tr, client_iris, netwk, statn,
                                loc_code, channel, sac_dict, t3, t4))
                else:
                    t_tr = threading.Thread(
                        target=lambda q, a, b, c, d, e, f, g, h:\
                            q.put(worker2(a, b, c, d, e, f, g, h)),
                        args = (cola_tr, client_gfz, netwk, statn,
                                loc_code, channel, sac_dict, t3, t4))
                t_mt = threading.Thread(
                    target=lambda q, a, b, c, d, e, f, g:\
                        q.put(worker4(a, b, c, d, e, f, g)),
                    args = (cola_mt, canal, netwk, statn, loc_code,
                            channel, t3, t4))
                t_tr.daemon = True
                t_tr.start()
                t_mt.daemon = True
                t_mt.start()
    time.sleep(100)
    return
    
    
def worker2(client, netwk, statn, loc_code, channel, sac_dict, time0, time1):
    try:
        st = client.get_waveforms(
            netwk, statn, loc_code, channel, time0, time1)
        name = '{}{}_{}{}.sac'.format(netwk, statn, channel, loc_code)
        if channel[:2] in ['HL', 'HN'] and len(st) > 1:
            return
        st[0].stats.sac = sac_dict 
        st.write(name, format='SAC', byteorder=0)
    except Exception as e:
        print(1, netwk, statn, loc_code, channel)
        print(e)
    return
    
    
def worker3(iris_client, netwk, statn, loc_code, channel, t1, t2):
    try:
        response_data = 'SAC_PZs_{}_{}_{}_{}'.format(
            netwk, statn, channel, loc_code) 
        if len(loc_code) is 0:
            response_data = 'SAC_PZs_{}_{}_{}___'.format(
                netwk, statn, channel)
        iris_client.sacpz(
            netwk, statn, loc_code, channel, t1, t2, filename=response_data)
    except:
        print(2, netwk, statn, loc_code, channel)
    return
        

def worker4(canal, netwk, statn, loc_code, channel, t1, t2):
    try:
        response = canal.response
        response_data = 'SAC_PZs_{}_{}_{}_{}'.format(
            netwk, statn, channel, loc_code) 
        if len(loc_code) is 0:
            response_data = 'SAC_PZs_{}_{}_{}___'.format(netwk, statn, channel)
        sacpz = response.get_sacpz()
        with open(response_data, 'w') as file:
            file.write('{}\n'.format(sacpz))
    except:
        print(3, netwk, statn, loc_code, channel)
    return


def __get_channel_information_manual(channel, lat_ep, lon_ep, depth):
    """Creamos un diccionario con la informacion para una estacion y una de las
    componentes, o canales, de dicha estacion.
    """
    station_lat = channel.latitude
    station_lon = channel.longitude
    loc_code = channel.location_code
    cmpaz = channel.azimuth
    cmpinc = 90 + channel.dip
    sac_dict = {
        'stla': station_lat,
        'stlo': station_lon,
        'khole': loc_code,
        'cmpaz': cmpaz,
        'cmpinc': cmpinc,
        'evla': lat_ep,
        'evlo': lon_ep,
        'evdp': depth
    }
    return sac_dict
    
    
if __name__ == '__main__':
    from obspy.core.utcdatetime import UTCDateTime
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-o", "--option",
                        choices=[
                            'acquisition',
                            'add_station',
                            'safe_acquisition',
                            'eew'
                        ], required=True, help="which method to run")
    parser.add_argument("-y", "--year", help="year of origin time")
    parser.add_argument("-m", "--month", help="month of origin time")
    parser.add_argument("-d", "--day", help="day of origin time")
    parser.add_argument("-ho", "--hour", help="hour of origin time")
    parser.add_argument("-mi", "--minute", help="minute of origin time")
    parser.add_argument("-s", "--second", help="second of origin time")
    parser.add_argument("--lat", help="latitude of hypocenter")
    parser.add_argument("--lon", help="longitude of hypocenter")
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.option == 'acquisition':
        event_time = UTCDateTime(2015, 9, 16, 22, 54, 33)
        lat_ep = -31
        lon_ep = -71
        depth = 25.0
        timedelta = 81 * 60
        time0 = time.time()
        acquisition(event_time, lat_ep, lon_ep, depth, ['strong'])
        print('time spent downloading metadata: ', time.time() - time0)
