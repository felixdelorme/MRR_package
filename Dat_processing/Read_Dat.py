import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
import os
from datetime import datetime

""" Félix DELORME on Tuesday 18 June 2024 """

def load(file):
    if not os.path.exists(file):
        raise FileNotFoundError(f"Le fichier {file} n'existe pas.")
    
    data = pd.read_csv(file, skiprows=[0, 2, 3], sep=',', parse_dates=[0], index_col=[0])
    data.index.name = 'Time'
    print("Noms des colonnes :", data.columns)
    print(data.head())
    
    return data

def verif_dates(data, event_start, event_stop):
    try:
        event_start = pd.to_datetime(event_start)
        event_stop = pd.to_datetime(event_stop)
    except ValueError:
        raise ValueError("Les dates fournies ne sont pas correctes.")
    
    if event_start < data.index.min() or event_stop > data.index.max():
        raise ValueError("Les dates fournies ne sont pas dans l'intervalle des données.")
    
    return data[event_start:event_stop]


def plot_temperature(data, filename):
    plt.figure(figsize=(10, 5))
    plt.plot_date(data.index, data['Prom_Temp_aire'], '-+', label='Average Air Temperature')
    plt.ylabel('Temperature (°C)')
    plt.xlabel('Time [UTC]')
    plt.title('Average Air Temperature - ' + filename)
    plt.grid()
    plt.legend()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%d/%m %H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.show()

def plot_precipitation_temperature(data_event, data_event_rain_hour, filename):
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.set_xlabel('Time [UTC]')
    ax1.set_ylabel('Precipitation intensity [mm/h]', color='blue')
    ax1.bar(data_event_rain_hour.index, data_event_rain_hour['Precip[mm/h]'], color='lightblue', edgecolor='blue', width=1/24)
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()
    ax2.set_ylabel('Temperature (°C)', color='green')
    ax2.plot(data_event.index, data_event['Prom_Temp_aire'], color='green', label='Average Air Temperature')
    ax2.tick_params(axis='y', labelcolor='green')

    plt.title('Precipitation and Temperature - ' + filename)
    fig.tight_layout()
    plt.grid(alpha=0.5)
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%d/%m %H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.show()

def plot_precipitation_albedo(data_event, data_event_rain_hour, filename):
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.set_xlabel('Time [UTC]')
    ax1.set_ylabel('Precipitation intensity [mm/h]', color='blue')
    ax1.bar(data_event_rain_hour.index, data_event_rain_hour['Precip[mm/h]'], color='lightblue', edgecolor='blue', width=1/24)
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()
    ax2.set_ylabel('Albedo', color='red')
    ax2.plot(data_event.index, data_event['Prom_rad_albedo'], color='red', label='Albedo')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylim(0, 1)

    plt.title('Precipitation and Albedo every 15min ' + filename)
    fig.tight_layout()
    plt.grid(alpha=0.5)
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%d/%m %H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.show()

def plot_daily_albedo(data):
    data_filtered = data[(data['Prom_rad_albedo'] >= 0) & (data['Prom_rad_albedo'] <= 1)]
    data_filtered['date'] = data_filtered.index.date
    daily_albedo = data_filtered.groupby('date')['Prom_rad_albedo'].mean().reset_index()
    daily_albedo.columns = ['Date', 'Albedo Journalier']
    daily_albedo.to_csv('daily_albedo.csv', index=False)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(daily_albedo['Date'], daily_albedo['Albedo Journalier'], marker='o', linestyle='-')
    ax.set_xlabel('Date')
    ax.set_ylabel('Albedo Journalier')
    ax.set_title('Albedo Journalier en fonction du temps')
    ax.grid(True)
    ax.set_ylim(0, 1)
    fig.autofmt_xdate()
    date_format = mdate.DateFormatter('%d/%m')
    ax.xaxis.set_major_formatter(date_format)
    plt.show()

####### add any function to view what you wish on a dat file ######

"""                               MAIN                                  """

def main():
    file = "/Users/felixdelorme/ESISAR/Stage/stage_felix_2024_MRR_Bolivie/2 Data/Plataforma Dat/2023Plataforma_.dat"
    filename = os.path.basename(file)
    data = load(file)

    event_start = '2023-11-01 00:00:00'
    event_stop = '2023-11-15 00:00:00'
    data_event = verif_dates(data, event_start, event_stop)
    
    data_event_rain_hour = data_event.resample(rule='1h').sum()
    data_event_rain_hour['Precip[mm/h]'] = data_event_rain_hour['Total_precipitacion']
    data_event_rain_hour['PrecipCum[mm]'] = data_event_rain_hour['Precip[mm/h]'].cumsum()
    data_event['Prom_Temp_aire'] = pd.to_numeric(data_event['Prom_Temp_aire'], errors='coerce')
    data_event = data_event.dropna(subset=['Prom_Temp_aire'])

    options = {
        "1": ("Afficher la température", plot_temperature),
        "2": ("Afficher les précipitations et la température", plot_precipitation_temperature),
        "3": ("Afficher les précipitations et l'albédo", plot_precipitation_albedo),
        "4": ("Afficher l'albédo journalier", plot_daily_albedo)
    }

    print("Options disponibles:")
    for key, (description, _) in options.items():
        print(f"{key}: {description}")

    choices = input("Choisissez une ou plusieurs options (séparées par des virgules): ").split(',')

    for choice in choices:
        choice = choice.strip()
        if choice in options:
            _, plot_function = options[choice]
            if choice == "1":
                plot_function(data_event, filename)
            elif choice == "2":
                plot_function(data_event, data_event_rain_hour, filename)
            elif choice == "3":
                plot_function(data_event, data_event_rain_hour, filename)
            elif choice == "4":
                plot_function(data_event)
        else:
            print(f"Option invalide: {choice}")

if __name__ == "__main__":
    main()
