import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import sys
from radiosonde.load_wyo_upperair import WyomingUpperAirSonde
from radiosonde.calc_wetbulb import wet_bulb_temperature
from plots.radiosonde.tephigram.tephigram_main import Tephigram


def plot_wyoming_tephigram(date, station, output_dir, theta_w=False):
    sonde_obj = WyomingUpperAirSonde(date, station)

    sonde_data = sonde_obj.get_dataframe()
    sonde_metadata = sonde_obj.get_metadata()
    pressures_winds = sonde_obj.prune_data(
        np.array([1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 300, 250, 200, 100, 50]))

    tpg = Tephigram()

    tpg.plot_profile(sonde_data['pressure'].values, sonde_data['temperature'].values,
                     label='Temperature', color='black', linewidth=1)
    tpg.plot_profile(sonde_data['pressure'].values, sonde_data['dewpoint'].values,
                     label='Dew Point', color='black', linewidth=1, linestyle='--')
    if theta_w:
        wet_bulb = wet_bulb_temperature(sonde_data['pressure'].values * 100,
                                        sonde_data['temperature'].values + 273.15,
                                        sonde_data['dewpoint'].values + 273.15)
        tpg.plot_profile(sonde_data['pressure'].values, wet_bulb - 273.15,
                         label='Wet Bulb', color='violet', linewidth=0.8)

    tpg.plot_barbs(pressures_winds['pressure'].values, pressures_winds['speed'].values,
                   pressures_winds['direction'].values + 180.0)
    tpg.plot_main_title(sonde_metadata)
    tpg.read_metadata(sonde_metadata)
    tpg.save_tephi(output_dir=output_dir)
    plt.close('all')


if __name__ == '__main__':

    # Change the following lines
    main_date = datetime(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
#    station_list = ['03005', '03023', '03238', '03354', '03502', '03590',
#                    '03693', '03743', '03808', '03882', '03918', '03953']
    station_list = [sys.argv[1]]
    output_dir = sys.argv[6]
    ############################

    for station_num in station_list:
        try:
            plot_wyoming_tephigram(date=main_date,
                                   station=station_num,
                                   output_dir=output_dir,
                                   theta_w=False)
            print(f"Plotted {main_date.strftime('%Y-%m-%d %HZ')} for station {station_num}.")
        except ValueError as ve:
            print(ve)
