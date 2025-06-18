#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 11:20:28 2025

@author: GeoSeisUtilities
"""

### This is the GUI for performing LET with FMTOMO-Extras

## Importing module
import glob
from obspy.core.event import read_events, Event, Origin, Magnitude
from obspy import Catalog
from obspy.core.event import Pick, Arrival, ResourceIdentifier
from obspy.core.utcdatetime import UTCDateTime
import pandas as pd
import pathlib
import numpy as np
import customtkinter as ctk
from tkinter import IntVar
import os

## Defining functions
def obspy2fmtomo(catalogue, stations, output_dir, phases):
    """
    Generate input files for the FMTOMO software package from an ObsPy
    catalogue.

    Parameters
    ----------
    catalogue : `obspy.Catalog` object
        Contains a list of `obspy.Event` objects, detailing origin times and
        phase picks.
    stations : `pandas.DataFrame` object
        DataFrame containing station information (latitude/longitude/elevation
        and a uid).
    output_dir : str
        Directory in which to save the output.
    phases : list of str
        List of phases to include in the output files.

    """

    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    pick_cols = ["olat", "olon", "odep", "ttime", "ttime_err"]
    pick_dict = {}
    for i, station in stations.iterrows():
        stat = station["Name"]
        print(f"Creating pick file for {stat}...")
        # Create DataFrames to store all picks for each phase
        dfs = {}
        for phase in phases:
            dfs[phase] = pd.DataFrame(columns=pick_cols)
        for event in catalogue:
            if event.preferred_origin() == None:
                origin = event.origins[0]
            else:
                origin = event.preferred_origin()
            olat, olon, odep = origin.latitude, origin.longitude, origin.depth
            otime = origin.time
            for pick in event.picks:
                if pick.waveform_id.station_code == stat:
                    line = pd.DataFrame([float(olat), float(olon), float(odep),
                                         float(pick.time - otime),
                                         float(pick.time_errors.uncertainty)],
                                        index=pick_cols).T
                    phase = pick.phase_hint
                    if phase in phases:
                        dfs[phase] = pd.concat([dfs[phase], line])
        for phase in phases:
            out = output_dir / "picks"
            out.mkdir(parents=True, exist_ok=True)
            outfile = out / f"{stat}.{phase}pick"
            out_df = dfs[phase]
            if out_df.empty:
                continue
            out_df = out_df.apply(pd.to_numeric)
            with outfile.open("w") as f:
                print(f"{len(out_df)}", file=f)
                for i, pick in out_df.iterrows():
                    print((f"{pick.olat:.4f} {pick.olon:.4f} "
                           f"{pick.odep/1000:.5f} {pick.ttime:.4f} "
                           f"{pick.ttime_err:.2f}"), file=f)
        anypicks = [dfs[phase].empty for pick in phases]
        if not np.all(anypicks):
            pick_dict[stat] = dfs
    print("Generation of pick files complete.")
    with (output_dir / "sourceslocal.in").open("w") as f:
        print(f"{len(pick_dict)}", file=f)
        for key, value in pick_dict.items():
            station = stations[stations["Name"] == key].iloc[0]
            stat = station["Name"]
            print((f"{station.Latitude:.4f} "
                   f"{station.Longitude:.4f} "
                   f"{station.Elevation:.4f}"), file=f)
            print(f"{len(phases)}", file=f)
            for phase in phases:
                print(f"1 1 {stat}.{phase}pick", file=f)
    stat_df = pd.DataFrame(columns=["Name", "Latitude", "Longitude", "Elevation"])
    for key in pick_dict.keys():
        stat_df = pd.concat([stat_df, stations[stations["Name"] == key]])
    stat_df.to_csv(str(output_dir / "stations_file.txt"), index=False)

class AlertWidget(ctk.CTkToplevel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.geometry("400x300")

        self.label = ctk.CTkLabel(self, text="ToplevelWindow")
        self.label.pack(padx=20, pady=20)

def convert_depth(depth_format, dataframe, column_name):
    if depth_format == "Elevation (m)":
        dataframe[column_name] = dataframe[column_name].apply(lambda x: -(float(x) / 1000))
    elif depth_format == "Elevation (Km)":
        dataframe[column_name] = dataframe[column_name].apply(lambda x: -float(x))
    elif depth_format == "Depth (m)":
        dataframe[column_name] = dataframe[column_name].apply(lambda x: float(x) / 1000)
    else:
        dataframe[column_name] = dataframe[column_name].apply(lambda x: float(x))

def convert_coord(coord_format, dataframe, column_name):
    if coord_format == "dd mm ss":
        dataframe[column_name] = dataframe[column_name].apply(lambda x: \
                                    int(x.split(' ')[0]) + (float(x.split(' ')[1])\
                                    +float(x.split(' ')[2])/60)/60)
    elif coord_format == "dd mm.mm":
        dataframe[column_name] = dataframe[column_name].apply(lambda x: \
                                    int(x.split(' ')[0]) + float(x.split(' ')[1])/60)
    elif coord_format == "ddLmm.mm":
        if 'N' in dataframe[column_name].iloc[0]:
            spac = 'N'
        elif 'n' in dataframe[column_name].iloc[0]:
            spac = 'n'
        elif 'S' in dataframe[column_name].iloc[0]:
            spac = 'S'
        elif 's' in dataframe[column_name].iloc[0]:
            spac = 's'
        elif 'E' in dataframe[column_name].iloc[0]:
            spac = 'E'
        elif 'e' in dataframe[column_name].iloc[0]:
            spac = 'e'
        elif 'W' in dataframe[column_name].iloc[0]:
            spac = 'W'
        elif 'w' in dataframe[column_name].iloc[0]:
            spac = 'w'
        else:
            return 'Wrong format for the longitude'
        dataframe[column_name] = dataframe[column_name].apply(lambda x: \
                                    int(x.split(spac)[0]) + float(x.split(spac)[1])/60)
    elif coord_format == "dd-mm-ss":
        dataframe[column_name] = dataframe[column_name].apply(lambda x: \
                                    int(x.split('-')[0]) + (float(x.split('-')[1])\
                                    +float(x.split('-')[2])/60)/60)
    elif coord_format == "dd L mm ss":
        if 'N' in dataframe[column_name].iloc[0]:
            spac = 'N'
        elif 'n' in dataframe[column_name].iloc[0]:
            spac = 'n'
        elif 'S' in dataframe[column_name].iloc[0]:
            spac = 'S'
        elif 's' in dataframe[column_name].iloc[0]:
            spac = 's'
        elif 'E' in dataframe[column_name].iloc[0]:
            spac = 'E'
        elif 'e' in dataframe[column_name].iloc[0]:
            spac = 'e'
        elif 'W' in dataframe[column_name].iloc[0]:
            spac = 'W'
        elif 'w' in dataframe[column_name].iloc[0]:
            spac = 'w'
        else:
            return 'Wrong format for the longitude'
        dataframe[column_name] = dataframe[column_name].apply(lambda x: \
                                    int(x.split(spac)[0]) + (float(x.split(spac)[2])\
                                    +float(x.split(spac)[3])/60)/60)
    elif coord_format == "ddLmmss.sss":
        if 'N' in dataframe[column_name].iloc[0]:
            spac = 'N'
        elif 'n' in dataframe[column_name].iloc[0]:
            spac = 'n'
        elif 'S' in dataframe[column_name].iloc[0]:
            spac = 'S'
        elif 's' in dataframe[column_name].iloc[0]:
            spac = 's'
        elif 'E' in dataframe[column_name].iloc[0]:
            spac = 'E'
        elif 'e' in dataframe[column_name].iloc[0]:
            spac = 'e'
        elif 'W' in dataframe[column_name].iloc[0]:
            spac = 'W'
        elif 'w' in dataframe[column_name].iloc[0]:
            spac = 'w'
        else:
            return 'Wrong format for the longitude'
        dataframe[column_name] = dataframe[column_name].apply(lambda x: \
                                    int(x.split(spac)[0]) + (float(x.split(spac)[1][0:-6])\
                                    +float(x.split(spac)[1][-6:])/60)/60)
    elif coord_format == "dd.dd":
        dataframe[column_name] = dataframe[column_name].apply(lambda x: float(x))
    else:
        return 'Wrong format for the longitude'

def parse_hypoellipse(file_path):
    picking_weight = {0:0.015, 1:0.03, 2:0.05, 3:0.08, 4:0.1}
    catalog = Catalog()
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("date") and "origin" in line:
            origin_line = lines[i + 1].strip()
            coord_line = lines[i + 2].strip()
            error_lines = lines[i + 4].strip()
            try:
                j = 8
                phases = []
                subsetline = lines[i + j].strip()
                phases.append(subsetline)
                while '-- magnitude data --' not in subsetline:
                    j += 1
                    subsetline = lines[i + j].strip()
                    phases.append(subsetline)
                phases = phases[:-2]
            except:
                print("No phases information")
            try:
                date_str = origin_line[:8]
                time_str = origin_line[9:19].replace(" ", "0")
                time = UTCDateTime(f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}T{time_str[:2]}:{time_str[2:4]}:{time_str[5:]}")
                lat = float(coord_line.split(' ')[0])
                lon = float(coord_line.split(' ')[-1])
                depth = float(origin_line[38:45].strip())
                try:
                    mag = float(origin_line[45:53].strip())
                    mag_ind = 0
                    magnitude = Magnitude(mag=mag)
                except:
                    mag_ind = 1
                if float(error_lines[:4].strip()) == 99.0:
                    seh = np.nan
                else:
                    seh = float(error_lines[:4].strip())
                sez = float(error_lines[4:9].strip())
                origin = Origin(time=time, latitude=lat, longitude=lon, depth=depth * 1000,
                        latitude_errors = seh, longitude_errors = seh, depth_errors = sez)
            except:
                print("Error in getting information from file")
            try:
                if mag_ind == 0:
                    event = Event(origins=[origin], magnitudes=[magnitude],
                                  resource_id=ResourceIdentifier(prefix="event"))
                else:
                    event = Event(origins=[origin],
                                  resource_id=ResourceIdentifier(prefix="event"))
                try:
                    dist_ = 0
                    azimuth_ = 0
                    for ph in phases:
                        station = ph[0:5].strip()
                        component = f"EH{ph[5].upper()}"
                        phase_hint = ph[12].strip().upper()
                        try:
                            time_weight = picking_weight[int(ph[14])]
                        except:
                            time_weight = picking_weight[4]
                        pick_time = float(ph[17:29].strip())
                        residual = float(ph[29:36].strip())
                        try:
                            std_error = float(ph[37:44].strip())
                        except:
                            std_error = np.nan
                        try:
                            dist = float(ph[44:51].strip())
                            dist_ = dist
                            azimuth = int(ph[51:56].strip())
                            azimuth_ = azimuth
                        except:
                            dist = dist_
                            azimuth = azimuth_
                        pick_time = time + pick_time
                        pick = Pick(time=pick_time, waveform_id={'station_code': station,
                            'channel_code': component}, phase_hint=phase_hint,
                            method_id="hypoellipse", time_errors=std_error)
                        arrival = Arrival(pick_id=pick.resource_id, phase=phase_hint,
                            azimuth=azimuth, distance=dist, time_residual=residual,
                            time_weight=time_weight)
                        event.picks.append(pick)
                        origin.arrivals.append(arrival)
                except:
                    print("Error in getting phases information")
                catalog.append(event)
            except:
                print("Error in adding event to catalog")
    return catalog

## Defining GUI
# general configuration
ctk.set_appearance_mode("System")
ctk.set_default_color_theme("green")
width = 800
height = 400
offset = 10

# main window
class GUI(ctk.CTk):
    def __init__(self):
        super().__init__()

        # configure window
        self.title("GUI LET FMTOMO")
        self.geometry(f"{width}x{height}")
        self.resizable(width=False, height=False)
        self.grid_rowconfigure((0, 1), weight=1)

        # general information
        self.entries = []
        self.buttons = []
        self.combobox = []
        self.fmtomo_folder = ctk.CTkEntry(self, placeholder_text="FMTOMO folder", width=width/4)
        self.fmtomo_folder.place(x=offset, y=offset)
        self.entries.append(self.fmtomo_folder)
        browse_button_1 = ctk.CTkButton(self, text="Browse", width=width/6,
                            command=lambda:[self.browse_folder(self.fmtomo_folder),
                                            self.change_color(browse_button_1)])
        browse_button_1.place(x=(2*offset + width/4), y=offset)
        self.buttons.append(browse_button_1)

        # general commands
        hitmap_button = ctk.CTkButton(self, text="Generate hitmap", width=width/4,
                            command=lambda:[self.perform_raytracing(),
                                            self.change_color(hitmap_button)])
        hitmap_button.place(x=width/2, y=offset)
        self.buttons.append(hitmap_button)
        reset_button = ctk.CTkButton(self, text="Reset GUI", width=width/7,
                            height=5*offset, command=lambda: self.reset_GUI(self.entries,
                                                                            self.buttons,
                                                                            self.combobox))
        reset_button.place(x=(offset + 5*(width/6)), y=offset)

        # create tabview
        self.tabview = ctk.CTkTabview(self, width=width, height=(height-(height/3)))
        self.tabview.grid(row=1, column=0, sticky="nsew")
        self.tabview.add("Earthquake data")
        self.tabview.add("Direct problem")
        self.tabview.add("Inverse problem")
        Eqs_data_tab = self.tabview.tab("Earthquake data")
        Dir_prob_tab = self.tabview.tab("Direct problem")
        Inv_prob_tab = self.tabview.tab("Inverse problem")

        # input selection tab
        eqs_switch = ctk.CTkSwitch(master=Eqs_data_tab, switch_height=10,
                            text="I need input generation",
                            command=lambda: self.switch_activation(eqs_widgets))
        eqs_switch.place(x=offset, y=0)
        self.cat_path = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Earthquakes catalog",
                                width=width/4)
        self.cat_path.place(x=width/2, y=6*offset)
        self.entries.append(self.cat_path)
        self.cat_path.configure(state="disabled")
        self.radio_var = IntVar(value=0)
        self.is_cat_folder = ctk.CTkRadioButton(Eqs_data_tab, text="Folder containing one file for earthquake",
                            variable=self.radio_var, value=0, radiobutton_height=9,
                            radiobutton_width=9, state="disabled")
        self.is_cat_folder.place(x=width/2, y=offset)
        self.is_cat_file = ctk.CTkRadioButton(Eqs_data_tab, text="File containing all catalog information",
                            variable=self.radio_var, value=1, radiobutton_height=9,
                            radiobutton_width=9, state="disabled")
        self.is_cat_file.place(x=width/2, y=3*offset)
        self.obspy_formats = {"QuakeML (qml)": "*.qml",
                         "Nlloc (hyp)": "*.hyp",
                         "SeisComP (xml)": "*.xml",
                         "Hypoellipse (.out)": "*.out",
                         "HypoDD (pha)": "*.pha",
                         "Seismic Handler (EVT)": "*.EVT",
                         "Seismic Handler (ASC)": "*.ASC",
                         "Seismic Handler (Q)": "*.Q",
                         "IASPEI Seismic Format (txt)": "*.txt",
                         "IASPEI Seismic Format (isf)": "*.isf",
                         "Nordic (out)": "*.out",
                         "ZMAP (txt)": "*.txt",
                         "Custom format": "*.*"}
        (("Text files","*.txt"), ("dat files","*.dat"),("csv files","*.csv"), ("All files","*.*"))
        
        self.browse_cat = ctk.CTkButton(Eqs_data_tab, text="Browse", width=width/6,
                            state="disabled", command=lambda: [self.browse_folder(self.cat_path)\
                            if self.radio_var.get() == 0 else self.browse_file(self.cat_path,\
                            list(self.obspy_formats.items())), self.change_color(self.browse_cat)])
        self.browse_cat.place(x=(offset + 3*width/4), y=6*offset)
        self.buttons.append(self.browse_cat)
        self.cat_format = ctk.CTkComboBox(Eqs_data_tab, width=width/4,
                                     values=list(self.obspy_formats.keys()),
                                     command=lambda x: self.second_switch_activator(\
                                     x, 'Custom format', custom_format_widgets))
        self.cat_format.set("Catalog format")
        self.combobox.append([self.cat_format, "Catalog format"])
        self.cat_format.place(x=width/2, y=9*offset)
        self.cat_format.configure(state="disabled")
        self.PS_choiche = ['P', 'S', 'P and S']
        self.PS = ctk.CTkComboBox(Eqs_data_tab, width=width/6,
                                     values=self.PS_choiche)
        self.PS.set("P or S")
        self.combobox.append([self.PS, "P or S"])
        self.PS.place(x=(3*width/4 + offset), y=9*offset)
        self.PS.configure(state="disabled")
        self.staz_path = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Station file",
                                width=width/4)
        self.staz_path.place(x=offset, y=3*offset)
        self.entries.append(self.staz_path)
        self.staz_path.configure(state="disabled")
        self.browse_staz = ctk.CTkButton(Eqs_data_tab, text="Browse", width=width/6,
                            state="disabled",  command=lambda:[self.browse_file(self.staz_path),
                                                               self.change_color(self.browse_staz)])
        self.browse_staz.place(x=(2*offset + width/4), y=3*offset)
        self.buttons.append(self.browse_staz)
        # station entries
        self.delm = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Delimiter",
                                width=width/4)
        self.delm.place(x=offset, y=6*offset)
        self.entries.append(self.delm)
        self.delm.configure(state="disabled")
        self.code_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Name (code) column",
                                width=width/4)
        self.code_col.place(x=offset, y=9*offset)
        self.entries.append(self.code_col)
        self.code_col.configure(state="disabled")
        self.coord_format = ["dd mm ss", "dd mm.mm", "ddLmm.mm", "dd-mm-ss",
                        "dd L mm ss", "ddLmmss.sss", "dd.dd"]
        self.lat_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Latitude column",
                                width=width/4)
        self.lat_col.place(x=offset, y=12*offset)
        self.entries.append(self.lat_col)
        self.lat_col.configure(state="disabled")
        self.lat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/6,
                                     values=self.coord_format)
        self.lat_form.set("Format")
        self.combobox.append([self.lat_form, "Format"])
        self.lat_form.place(x=(width/4 + 2*offset), y=12*offset)
        self.lat_form.configure(state="disabled")
        self.lon_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Longitude column",
                                width=width/4)
        self.lon_col.place(x=offset, y=15*offset)
        self.entries.append(self.lon_col)
        self.lon_col.configure(state="disabled")
        self.lon_form = ctk.CTkComboBox(Eqs_data_tab, width=width/6,
                                     values=self.coord_format)
        self.lon_form.set("Format")
        self.combobox.append([self.lon_form, "Format"])
        self.lon_form.place(x=(width/4 + 2*offset), y=15*offset)
        self.lon_form.configure(state="disabled")
        self.depth_format = ["Elevation (m)", "Depth (m)", "Elevation (Km)", "Depth (Km)"]
        self.elev_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Elevation column",
                                width=width/4)
        self.elev_col.place(x=offset, y=18*offset)
        self.entries.append(self.elev_col)
        self.elev_col.configure(state="disabled")
        self.elev_form = ctk.CTkComboBox(Eqs_data_tab, width=width/6,
                                     values=self.depth_format)
        self.elev_form.set("Format")
        self.combobox.append([self.elev_form, "Format"])
        self.elev_form.place(x=(width/4 + 2*offset), y=18*offset)
        self.elev_form.configure(state="disabled")
        # catalog entries
        self.delm_cat = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Cat. Delim.",
                                width=(width/8 - offset))
        self.delm_cat.place(x=width/2, y=13*offset)
        self.entries.append(self.delm_cat)
        self.delm_cat.configure(state="disabled")
        date_format={
                    "yyyy-mm-dd hh:mm:ss": "%Y-%m-%d %H:%M:%S",
                    "yyyy/mm/dd hh:mm:ss": "%Y/%m/%d %H:%M:%S",
                    "dd-mm-yyyy hh:mm:ss": "%d-%m-%Y %H:%M:%S",
                    "dd/mm/yyyy hh:mm:ss": "%d/%m/%Y %H:%M:%S",
                    "yyyy-mm-dd": "%Y-%m-%d",
                    "yyyy/mm/dd": "%Y/%m/%d",
                    "dd-mm-yyyy": "%d-%m-%Y",
                    "dd/mm/yyyy": "%d/%m/%Y",
                    "yyyy-mm-ddThh:mm:ss": "%Y-%m-%dT%H:%M:%S",
                    "yyyy-mm-ddThh:mm:ss.ss": "%Y-%m-%dT%H:%M:%S.%f",
                    "yyyy-mm-ddThh:mm:ssZ": "%Y-%m-%dT%H:%M:%SZ",
                    "yyyy-mm-ddThh:mm:ss.sssZ": "%Y-%m-%dT%H:%M:%S.%fZ",
                    "yyyymmddhhmmss": "%Y%m%d%H%M%S",
                    "yyyymmdd": "%Y%m%d",
                    "dd month yyyy hh:mm:ss": "%d %B %Y %H:%M:%S",
                    "day, dd month yyyy": "%A, %d %B %Y"
                    }
        self.time_cat_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Time col.",
                                width=(width/8 - offset))
        self.time_cat_col.place(x=width/2, y=16*offset)
        self.entries.append(self.time_cat_col)
        self.time_cat_col.configure(state="disabled")
        self.time_cat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/8,
                                     values=list(date_format.keys()))
        self.time_cat_form.set("Format")
        self.combobox.append([self.time_cat_form, "Format"])
        self.time_cat_form.place(x=(5*width/8), y=16*offset)
        self.time_cat_form.configure(state="disabled")
        self.lat_cat_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Lat. col.",
                                width=(width/8 - offset))
        self.lat_cat_col.place(x=width/2, y=19*offset)
        self.entries.append(self.lat_cat_col)
        self.lat_cat_col.configure(state="disabled")
        self.lat_cat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/8,
                                     values=self.coord_format)
        self.lat_cat_form.set("Format")
        self.combobox.append([self.lat_cat_form, "Format"])
        self.lat_cat_form.place(x=(5*width/8), y=19*offset)
        self.lat_cat_form.configure(state="disabled")
        self.lon_cat_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Lon. col.",
                                width=(width/8 - offset))
        self.lon_cat_col.place(x=width/2, y=22*offset)
        self.entries.append(self.lon_cat_col)
        self.lon_cat_col.configure(state="disabled")
        self.lon_cat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/8,
                                     values=self.coord_format)
        self.lon_cat_form.set("Format")
        self.combobox.append([self.lon_cat_form, "Format"])
        self.lon_cat_form.place(x=(5*width/8), y=22*offset)
        self.lon_cat_form.configure(state="disabled")
        self.elev_cat_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Dep. col.",
                                width=(width/8 - offset))
        self.elev_cat_col.place(x=width/2, y=25*offset)
        self.entries.append(self.elev_cat_col)
        self.elev_cat_col.configure(state="disabled")
        self.elev_cat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/8,
                                     values=self.depth_format)
        self.elev_cat_form.set("Format")
        self.combobox.append([self.elev_cat_form, "Format"])
        self.elev_cat_form.place(x=(5*width/8), y=25*offset)
        self.elev_cat_form.configure(state="disabled")
        # picks entries
        self.delm_pick = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Pick Delim.",
                                width=(width/8 - offset))
        self.delm_pick.place(x=(3*width/4 + offset), y=13*offset)
        self.entries.append(self.delm_pick)
        self.delm_pick.configure(state="disabled")
        picking_format = ["Still working on...", "Still working on..."]
        self.pick_form = ctk.CTkComboBox(Eqs_data_tab, width=width/9,
                                     values=picking_format)
        self.pick_form.set("Format")
        self.combobox.append([self.pick_form, "Format"])
        self.pick_form.place(x=(3*width/4 + offset), y=16*offset)
        self.pick_form.configure(state="disabled")
        # button
        self.gen_input = ctk.CTkButton(Eqs_data_tab, text="Generate input files", width=(2*(width/5) + 2*offset),
                             state='disabled', command=lambda:[self.generate_input(),
                                                               self.change_color(self.gen_input)])
        self.gen_input.place(x=offset, y=25*offset)
        self.buttons.append(self.gen_input)
        eqs_widgets = [self.cat_path, self.is_cat_folder, self.is_cat_file, self.cat_format, self.staz_path,
                       self.browse_cat, self.browse_staz, self.delm, self.code_col, self.lat_col, self.lat_form,
                       self.lon_col, self.lon_form, self.elev_col, self.elev_form, self.PS, self.gen_input]
        custom_format_widgets = []
        # custom_format_widgets = [self.delm_cat, self.time_cat_col, self.time_cat_form,
        #             self.lat_cat_col, self.lat_cat_form, self.lon_cat_col, self.lon_cat_form,
        #             self.elev_cat_col, self.elev_cat_form, self.delm_pick, self.pick_form]

        # direct problem tab
        dir_switch = ctk.CTkSwitch(master=Dir_prob_tab, switch_height=10,
                            text="I need to define the model",
                            command=lambda: self.switch_activation(dir_widgets))
        dir_switch.place(x=offset, y=0)
        # model entries
        self.min_lat = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Minimum Latitude",
                                width=(width/5 - 2*offset))
        self.min_lat.place(x=offset, y=3*offset)
        self.entries.append(self.min_lat)
        self.min_lat.configure(state="disabled")
        self.max_lat = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Maximum Latitude",
                                width=(width/5 - 2*offset))
        self.max_lat.place(x=(width/5), y=3*offset)
        self.entries.append(self.max_lat)
        self.max_lat.configure(state="disabled")
        P_lab = ctk.CTkLabel(Dir_prob_tab, text="Nodes P model")
        P_lab.place(x=(2*(width/5) - offset), y=0)
        S_lab = ctk.CTkLabel(Dir_prob_tab, text="Nodes S model")
        S_lab.place(x=(3*(width/5) - offset), y=0)
        self.n_points_latP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_latP.place(x=(2*(width/5) - offset), y=3*offset)
        self.entries.append(self.n_points_latP)
        self.n_points_latP.configure(state="disabled")
        self.n2_points_latP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_latP.place(x=(2.5*(width/5) - offset), y=3*offset)
        self.entries.append(self.n2_points_latP)
        self.n2_points_latP.configure(state="disabled")
        self.n_points_latS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_latS.place(x=(3*(width/5) - offset), y=3*offset)
        self.entries.append(self.n_points_latS)
        self.n_points_latS.configure(state="disabled")
        self.n2_points_latS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_latS.place(x=(3.5*(width/5) - offset), y=3*offset)
        self.entries.append(self.n2_points_latS)
        self.n2_points_latS.configure(state="disabled")
        self.min_lon = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Min. Longitude",
                                width=(width/5 - 2*offset))
        self.min_lon.place(x=offset, y=6*offset)
        self.entries.append(self.min_lon)
        self.min_lon.configure(state="disabled")
        self.max_lon = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Max. Longitude",
                                width=(width/5 - 2*offset))
        self.max_lon.place(x=(width/5), y=6*offset)
        self.entries.append(self.max_lon)
        self.max_lon.configure(state="disabled")
        self.n_points_lonP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_lonP.place(x=(2*(width/5) - offset), y=6*offset)
        self.entries.append(self.n_points_lonP)
        self.n_points_lonP.configure(state="disabled")
        self.n2_points_lonP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_lonP.place(x=(2.5*(width/5) - offset), y=6*offset)
        self.entries.append(self.n2_points_lonP)
        self.n2_points_lonP.configure(state="disabled")
        self.n_points_lonS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_lonS.place(x=(3*(width/5) - offset), y=6*offset)
        self.entries.append(self.n_points_lonS)
        self.n_points_lonS.configure(state="disabled")
        self.n2_points_lonS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_lonS.place(x=(3.5*(width/5) - offset), y=6*offset)
        self.entries.append(self.n2_points_lonS)
        self.n2_points_lonS.configure(state="disabled")
        self.min_dep = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Min. Elevation (<0)",
                                width=(width/5 - 2*offset))
        self.min_dep.place(x=offset, y=9*offset)
        self.entries.append(self.min_dep)
        self.min_dep.configure(state="disabled")
        self.max_dep = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Max. Depth (>0)",
                                width=(width/5 - 2*offset))
        self.max_dep.place(x=(width/5), y=9*offset)
        self.entries.append(self.max_dep)
        self.max_dep.configure(state="disabled")
        self.n_points_depP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_depP.place(x=(2*(width/5) - offset), y=9*offset)
        self.entries.append(self.n_points_depP)
        self.n_points_depP.configure(state="disabled")
        self.n2_points_depP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_depP.place(x=(2.5*(width/5) - offset), y=9*offset)
        self.entries.append(self.n2_points_depP)
        self.n2_points_depP.configure(state="disabled")
        self.n_points_depS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_depS.place(x=(3*(width/5) - offset), y=9*offset)
        self.entries.append(self.n_points_depS)
        self.n_points_depS.configure(state="disabled")
        self.n2_points_depS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_depS.place(x=(3.5*(width/5) - offset), y=9*offset)
        self.entries.append(self.n2_points_depS)
        self.n2_points_depS.configure(state="disabled")
        self.refined_grid = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Refine grid dim.",
                                width=(width/5 - 3*offset))
        self.refined_grid.configure(state="disabled")
        self.entries.append(self.refined_grid)
        self.refined_grid.place(x=(4*(width/5) - offset), y=4.5*offset)
        self.n_cell_refined = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Refine grid ext.",
                                width=(width/5 - 3*offset))
        self.n_cell_refined.configure(state="disabled")
        self.entries.append(self.n_cell_refined)
        self.n_cell_refined.place(x=(4*(width/5) - offset), y=7.5*offset)
        # estimate spacing
        self.est_dist = ctk.CTkButton(Dir_prob_tab, text="Estimate point spacing", width=(width - 5*offset),
                             state='disabled', command=lambda:[self.estimate_spacing(),
                                                               self.change_color(self.est_dist)])
        self.est_dist.place(x=offset, y=13*offset)
        self.buttons.append(self.est_dist)
        estimated_prop = ctk.CTkLabel(Dir_prob_tab, text="Propagation grid spacing (km)")
        estimated_prop.place(x=offset, y=16*offset)
        self.est_latP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lat. spac.",
                                width=(width/6 - 2*offset))
        self.est_latP.place(x=offset, y=19*offset)
        self.entries.append(self.est_latP)
        self.est_latP.configure(state="disabled")
        self.est_lonP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lon. spac.",
                                width=(width/6 - 2*offset))
        self.est_lonP.place(x=offset, y=22*offset)
        self.entries.append(self.est_lonP)
        self.est_lonP.configure(state="disabled")
        self.est_depP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Dep. spac.",
                                width=(width/6 - 2*offset))
        self.est_depP.place(x=offset, y=25*offset)
        self.entries.append(self.est_depP)
        self.est_depP.configure(state="disabled")
        self.est_latS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lat. spac.",
                                width=(width/6 - 2*offset))
        self.est_latS.place(x=width/6, y=19*offset)
        self.entries.append(self.est_latS)
        self.est_latS.configure(state="disabled")
        self.est_lonS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lon. spac.",
                                width=(width/6 - 2*offset))
        self.est_lonS.place(x=width/6, y=22*offset)
        self.entries.append(self.est_lonS)
        self.est_lonS.configure(state="disabled")
        self.est_depS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Dep. spac.",
                                width=(width/6 - 2*offset))
        self.est_depS.place(x=width/6, y=25*offset)
        self.entries.append(self.est_depS)
        self.est_depS.configure(state="disabled")
        estimated_vel = ctk.CTkLabel(Dir_prob_tab, text="Velocity grid spacing (km)")
        estimated_vel.place(x=width/3, y=16*offset)
        self.est_lat2P = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lat. spac.",
                                width=(width/6 - 2*offset))
        self.est_lat2P.place(x=width/3, y=19*offset)
        self.entries.append(self.est_lat2P)
        self.est_lat2P.configure(state="disabled")
        self.est_lon2P = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lon. spac.",
                                width=(width/6 - 2*offset))
        self.est_lon2P.place(x=width/3, y=22*offset)
        self.entries.append(self.est_lon2P)
        self.est_lon2P.configure(state="disabled")
        self.est_dep2P = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Dep. spac.",
                                width=(width/6 - 2*offset))
        self.est_dep2P.place(x=width/3, y=25*offset)
        self.entries.append(self.est_dep2P)
        self.est_dep2P.configure(state="disabled")
        self.est_lat2S = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lat. spac.",
                                width=(width/6 - 2*offset))
        self.est_lat2S.place(x=(width/2 - offset), y=19*offset)
        self.entries.append(self.est_lat2S)
        self.est_lat2S.configure(state="disabled")
        self.est_lon2S = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lon. spac.",
                                width=(width/6 - 2*offset))
        self.est_lon2S.place(x=(width/2 - offset), y=22*offset)
        self.entries.append(self.est_lon2S)
        self.est_lon2S.configure(state="disabled")
        self.est_dep2S = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Dep. spac.",
                                width=(width/6 - 2*offset))
        self.est_dep2S.place(x=(width/2 - offset), y=25*offset)
        self.entries.append(self.est_dep2S)
        self.est_dep2S.configure(state="disabled")
        tot_n_points = ctk.CTkLabel(Dir_prob_tab, text="Number of propagation nodes")
        tot_n_points.place(x=(2*(width/3) - offset), y=16*offset)
        self.est_tot_pointsP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="P nodes",
                                width=(width/6 - 2*offset))
        self.est_tot_pointsP.configure(state="disabled")
        self.entries.append(self.est_tot_pointsP)
        self.est_tot_pointsP.place(x=(2*(width/3) - offset), y=19*offset)
        self.est_tot_pointsS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="S nodes",
                                width=(width/6 - 2*offset))
        self.est_tot_pointsS.configure(state="disabled")
        self.entries.append(self.est_tot_pointsS)
        self.est_tot_pointsS.place(x=(5*(width/6) - 2*offset), y=19*offset)
        tot_n_points2 = ctk.CTkLabel(Dir_prob_tab, text="Number of velocity nodes")
        tot_n_points2.place(x=2*(width/3) - offset, y=22*offset)
        self.est_tot_points2P = ctk.CTkEntry(Dir_prob_tab, placeholder_text="P nodes",
                                width=(width/6 - 2*offset))
        self.est_tot_points2P.configure(state="disabled")
        self.entries.append(self.est_tot_points2P)
        self.est_tot_points2P.place(x=(2*(width/3) - offset), y=25*offset)
        self.est_tot_points2S = ctk.CTkEntry(Dir_prob_tab, placeholder_text="S nodes",
                                width=(width/6 - 2*offset))
        self.est_tot_points2S.configure(state="disabled")
        self.entries.append(self.est_tot_points2S)
        self.est_tot_points2S.place(x=(5*(width/6) - 2*offset), y=25*offset)
        dir_widgets = [self.min_lat, self.max_lat, self.n_points_latP, self.min_lon,
                       self.max_lon, self.n_points_lonP, self.min_dep, self.max_dep,
                       self.n_points_depP, self.n2_points_latP, self.n2_points_lonP,
                       self.n2_points_depP, self.n_points_latS, self.n_points_lonS,
                       self.n_points_depS, self.n2_points_latS, self.n2_points_lonS,
                       self.n2_points_depS,self.refined_grid, self.n_cell_refined,
                       self.est_dist]

        # inverse problem tab
        inv_switch = ctk.CTkSwitch(master=Inv_prob_tab, switch_height=10,
                            text="I need to set the inversion",
                            command=lambda: self.switch_activation(inv_widgets))
        inv_switch.place(x=offset, y=0)
        # inversion parameters
        self.n_iter = ctk.CTkEntry(Inv_prob_tab, placeholder_text="N. of iteration",
                                width=(width/5))
        self.n_iter.configure(state="disabled")
        self.entries.append(self.n_iter)
        self.n_iter.place(x=offset, y=3*offset)
        self.n_subs = ctk.CTkEntry(Inv_prob_tab, placeholder_text="N. of subspaces",
                                width=(width/5))
        self.n_subs.configure(state="disabled")
        self.entries.append(self.n_subs)
        self.n_subs.place(x=(width/5 + 2*offset), y=3*offset)
        self.dampP = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Damping value P",
                                width=(width/5))
        self.dampP.configure(state="disabled")
        self.entries.append(self.dampP)
        self.dampP.place(x=offset, y=6*offset)
        self.smoothP = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Smoothing value P",
                                width=(width/5))
        self.smoothP.configure(state="disabled")
        self.entries.append(self.smoothP)
        self.smoothP.place(x=offset, y=9*offset)
        self.dampS = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Damping value S",
                                width=(width/5))
        self.dampS.configure(state="disabled")
        self.entries.append(self.dampS)
        self.dampS.place(x=(width/5 + 2*offset), y=6*offset)
        self.smoothS = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Smoothing value S",
                                width=(width/5))
        self.smoothS.configure(state="disabled")
        self.entries.append(self.smoothS)
        self.smoothS.place(x=(width/5 + 2*offset), y=9*offset)
        # velocity model
        self.radio_var2 = IntVar(value=0)
        self.is_vel_file = ctk.CTkRadioButton(Inv_prob_tab, text="Velocity model is in a file",
                            variable=self.radio_var2, value=0, radiobutton_height=9,
                            radiobutton_width=9, state="disabled", command=lambda: \
                            self.third_switch_activator([self.vel_file, self.browse_vel],
                            [self.min_velP, self.min_velS, self.max_velP, self.max_velS]))
        self.is_vel_file.place(x=width/2, y=offset)
        self.is_vel_gradual = ctk.CTkRadioButton(Inv_prob_tab, text="Velocity is a gradual model",
                            variable=self.radio_var2, value=1, radiobutton_height=9,
                            radiobutton_width=9, state="disabled", command=lambda:\
                            self.third_switch_activator([self.min_velP, self.min_velS,
                            self.max_velP, self.max_velS], [self.vel_file, self.browse_vel]))
        self.is_vel_gradual.place(x=width/2, y=3*offset)
        self.vel_file = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Velocity file",
                                width=(width/4+offset))
        self.vel_file.place(x=width/2, y=6*offset)
        self.entries.append(self.vel_file)
        self.vel_file.configure(state="disabled")
        self.vel_filetypes = (("Velocity files","*.vel"),("Text files","*.txt"),
                              ("Dat files","*.dat"),("csv files","*.csv"),("All files","*.*"))
        self.browse_vel = ctk.CTkButton(Inv_prob_tab, text="Browse", width=width/6,
                            state="disabled",  command=lambda:[self.browse_file(self.vel_file,
                            self.vel_filetypes), self.change_color(self.browse_vel)])
        self.browse_vel.place(x=(2*offset + 3*(width/4)), y=6*offset)
        self.buttons.append(self.browse_vel)
        self.min_velP = ctk.CTkEntry(Inv_prob_tab, placeholder_text="P velocity top",
                                width=(width/5+offset))
        self.min_velP.place(x=width/2, y=10*offset)
        self.entries.append(self.min_velP)
        self.min_velP.configure(state="disabled")
        self.min_velS = ctk.CTkEntry(Inv_prob_tab, placeholder_text="S velocity top",
                                width=(width/5+offset))
        self.min_velS.place(x=(2*offset + 7*(width/10)), y=10*offset)
        self.entries.append(self.min_velS)
        self.min_velS.configure(state="disabled")
        self.max_velP = ctk.CTkEntry(Inv_prob_tab, placeholder_text="P velocity bottom",
                                width=(width/5+offset))
        self.max_velP.place(x=width/2, y=13*offset)
        self.entries.append(self.max_velP)
        self.max_velP.configure(state="disabled")
        self.max_velS = ctk.CTkEntry(Inv_prob_tab, placeholder_text="S velocity bottom",
                                width=(width/5+offset))
        self.max_velS.place(x=(2*offset + 7*(width/10)), y=13*offset)
        self.entries.append(self.max_velS)
        self.max_velS.configure(state="disabled")
        # synthetic models
        self.spike_switch = ctk.CTkSwitch(master=Inv_prob_tab, switch_height=10,
                            text="Build spike model", state="disabled",
                            command=lambda: self.third_switch_activator([self.pert_value,
                            self.spike_lat, self.spike_lon, self.spike_dep], [self.pert_value2,
                            self.check_size]) if self.spike_switch.get() == 1 else
                            self.third_switch_activator([self.pert_value2, self.check_size],
                            [self.pert_value, self.spike_lat, self.spike_lon, self.spike_dep]))
        self.spike_switch.place(x=offset, y=18*offset)
        self.pert_value = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Spike value",
                                width=(width/5))
        self.pert_value.place(x=offset, y=21*offset)
        self.entries.append(self.pert_value)
        self.pert_value.configure(state="disabled")
        self.spike_lat = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Lat",
                                width=(width/6 - 3*offset))
        self.spike_lat.place(x=offset, y=24*offset)
        self.entries.append(self.spike_lat)
        self.spike_lat.configure(state="disabled")
        self.spike_lon = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Lon",
                                width=(width/6 - 3*offset))
        self.spike_lon.place(x=(width/6 - offset), y=24*offset)
        self.entries.append(self.spike_lon)
        self.spike_lon.configure(state="disabled")
        self.spike_dep = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Dep",
                                width=(width/6 - 3*offset))
        self.spike_dep.place(x=(width/3 - 3*offset), y=24*offset)
        self.entries.append(self.spike_dep)
        self.spike_dep.configure(state="disabled")
        self.check_switch = ctk.CTkSwitch(master=Inv_prob_tab, switch_height=10,
                            text="Build checkerboard model", state="disabled",
                            command=lambda: self.third_switch_activator([self.pert_value2,
                            self.check_size], [self.pert_value, self.spike_lat,
                            self.spike_lon, self.spike_dep]) if self.check_switch.get() == 1
                            else self.third_switch_activator([self.pert_value,
                            self.spike_lat, self.spike_lon, self.spike_dep], [self.pert_value2,
                            self.check_size]))
        self.check_switch.place(x=width/2, y=18*offset)
        self.pert_value2 = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Perturbation value",
                                width=(width/5+offset))
        self.pert_value2.place(x=width/2, y=21*offset)
        self.entries.append(self.pert_value2)
        self.pert_value2.configure(state="disabled")
        self.check_size = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Check. dimension",
                                width=(width/5+offset))
        self.check_size.place(x=width/2, y=24*offset)
        self.entries.append(self.check_size)
        self.check_size.configure(state="disabled")
        self.comp_inp = ctk.CTkButton(Inv_prob_tab, text="Compile FMTOMO\n files",
                            width=(width/5 + offset), height=(height/6 - offset),
                            state='disabled', command=lambda:[self.compile_input_files(),
                                                              self.change_color(self.comp_inp)])
        self.comp_inp.place(x=(2*offset + 7*(width/10)), y=21*offset)
        self.buttons.append(self.comp_inp)
        inv_widgets = [self.dampP, self.smoothP, self.dampS, self.smoothS, self.n_iter,
                       self.is_vel_file, self.is_vel_gradual, self.vel_file, self.comp_inp,
                       self.n_subs, self.browse_vel, self.spike_switch, self.check_switch]

    def browse_folder(self, entry):
        from tkinter import filedialog
        filename = filedialog.askdirectory(initialdir = "/home/", title = "Browse folders",)
        entry.delete(0,)
        entry.insert(0,str(filename))

    def browse_file(self, entry, filetypes=(("Text files","*.txt"),
                                            ("dat files","*.dat"),("csv files","*.csv"),
                                            ("All files","*.*"))):
        from tkinter import filedialog
        filename=filedialog.askopenfilename(initialdir = "/home/",title = "Open file",
                                            filetypes = filetypes)
        entry.delete(0,)
        entry.insert(0,str(filename))

    def change_color(self, button):
        button.configure(fg_color="blue", hover_color="darkblue")

    def switch_activation(self, widgets):
        state = widgets[0].cget("state")
        if state == "normal":
            state = "disabled"
        else:
            state = "normal"
        for wid in widgets:
            wid.configure(state=state)

    def second_switch_activator(self, entry, value, widgets):
        if entry == value:
            self.switch_activation(widgets)
        else:
            try:
                state = widgets[0].cget("state")
                if state != "disabled":
                    self.switch_activation(widgets)
            except:
                pass

    def third_switch_activator(self, activate, deactivate):
        for wid in activate:
            wid.configure(state="normal")
        for wid2 in deactivate:
            wid2.configure(state="disabled")

    def generate_input(self):
        self.checked_eqs = 1
        self.resize_staz_eqs()
        # acquire catalog info
        fmtomo_path = self.fmtomo_folder.get()
        ps = self.PS.get()
        if ps == 'P':
            obspy2fmtomo(self.ev_cat, self.staz_df, f'{fmtomo_path}/invert_p/mkdir/', ["P"])
        elif ps == 'S':
            obspy2fmtomo(self.ev_cat, self.staz_df, f'{fmtomo_path}/invert_s/mkdir/', ["S"])
        elif ps == 'P and S':
            obspy2fmtomo(self.ev_cat, self.staz_df, f'{fmtomo_path}/invert_p/mkdir/', ["P"])
            obspy2fmtomo(self.ev_cat, self.staz_df, f'{fmtomo_path}/invert_s/mkdir/', ["S"])

    def estimate_spacing(self):
        from geopy.distance import geodesic
        widgets = [self.est_latP, self.est_lat2P, self.est_lonP, self.est_lon2P,
                   self.est_depP, self.est_dep2P, self.est_latS, self.est_lat2S,
                   self.est_lonS, self.est_lon2S, self.est_depS, self.est_dep2S,
                   self.est_tot_pointsP, self.est_tot_points2P,
                   self.est_tot_pointsS, self.est_tot_points2S]
        self.switch_activation(widgets)
        lat1, lon1 = float(self.min_lat.get()), float(self.min_lon.get())
        lat2, lon2 = float(self.max_lat.get()), float(self.max_lon.get())
        dep1, dep2 = float(self.min_dep.get()), float(self.max_dep.get())
        nlatP, nlonP, ndepP = int(self.n_points_latP.get())+1, int(self.n_points_lonP.get())+1, int(self.n_points_depP.get())+1
        nlatS, nlonS, ndepS = int(self.n_points_latS.get())+1, int(self.n_points_lonS.get())+1, int(self.n_points_depS.get())+1
        nlat2P, nlon2P, ndep2P = int(self.n2_points_latP.get())+1, int(self.n2_points_lonP.get())+1, int(self.n2_points_depP.get())+1
        nlat2S, nlon2S, ndep2S = int(self.n2_points_latS.get())+1, int(self.n2_points_lonS.get())+1, int(self.n2_points_depS.get())+1
        lat_spacP = geodesic((lat1, lon1), (lat2, lon1)).kilometers / nlatP
        lon_spacP = geodesic((lat1, lon1), (lat1, lon2)).kilometers / nlonP
        dep_spacP = abs(dep2 - dep1) / ndepP
        lat_spac2P = geodesic((lat1, lon1), (lat2, lon1)).kilometers / nlat2P
        lon_spac2P = geodesic((lat1, lon1), (lat1, lon2)).kilometers / nlon2P
        dep_spac2P = abs(dep2 - dep1) / ndep2P
        lat_spacS = geodesic((lat1, lon1), (lat2, lon1)).kilometers / nlatS
        lon_spacS = geodesic((lat1, lon1), (lat1, lon2)).kilometers / nlonS
        dep_spacS = abs(dep2 - dep1) / ndepS
        lat_spac2S = geodesic((lat1, lon1), (lat2, lon1)).kilometers / nlat2S
        lon_spac2S = geodesic((lat1, lon1), (lat1, lon2)).kilometers / nlon2S
        dep_spac2S = abs(dep2 - dep1) / ndep2S
        while self.est_latP.get() != '':
            self.est_latP.delete(0,)
        self.est_latP.insert(0,str(round(lat_spacP,3)))
        while self.est_lonP.get() != '':
            self.est_lonP.delete(0,)
        self.est_lonP.insert(0,str(round(lon_spacP,3)))
        while self.est_depP.get() != '':
            self.est_depP.delete(0,)
        self.est_depP.insert(0,str(round(dep_spacP,3)))
        while self.est_lat2P.get() != '':
            self.est_lat2P.delete(0,)
        self.est_lat2P.insert(0,str(round(lat_spac2P,3)))
        while self.est_lon2P.get() != '':
            self.est_lon2P.delete(0,)
        self.est_lon2P.insert(0,str(round(lon_spac2P,3)))
        while self.est_dep2P.get() != '':
            self.est_dep2P.delete(0,)
        self.est_dep2P.insert(0,str(round(dep_spac2P,3)))
        while self.est_latS.get() != '':
            self.est_latS.delete(0,)
        self.est_latS.insert(0,str(round(lat_spacS,3)))
        while self.est_lonS.get() != '':
            self.est_lonS.delete(0,)
        self.est_lonS.insert(0,str(round(lon_spacS,3)))
        while self.est_depS.get() != '':
            self.est_depS.delete(0,)
        self.est_depS.insert(0,str(round(dep_spacS,3)))
        while self.est_lat2S.get() != '':
            self.est_lat2S.delete(0,)
        self.est_lat2S.insert(0,str(round(lat_spac2S,3)))
        while self.est_lon2S.get() != '':
            self.est_lon2S.delete(0,)
        self.est_lon2S.insert(0,str(round(lon_spac2S,3)))
        while self.est_dep2S.get() != '':
            self.est_dep2S.delete(0,)
        self.est_dep2S.insert(0,str(round(dep_spac2S,3)))
        while self.est_tot_pointsP.get() != '':
            self.est_tot_pointsP.delete(0,)
        self.est_tot_pointsP.insert(0,str(nlatP*nlonP*ndepP))
        while self.est_tot_points2P.get() != '':
            self.est_tot_points2P.delete(0,)
        self.est_tot_points2P.insert(0,str(nlat2P*nlon2P*ndep2P))
        while self.est_tot_pointsS.get() != '':
            self.est_tot_pointsS.delete(0,)
        self.est_tot_pointsS.insert(0,str(nlatS*nlonS*ndepS))
        while self.est_tot_points2S.get() != '':
            self.est_tot_points2S.delete(0,)
        self.est_tot_points2S.insert(0,str(nlat2S*nlon2S*ndep2S))
        self.switch_activation(widgets)

    def resize_staz_eqs(self):
        lat1, lon1 = float(self.min_lat.get()), float(self.min_lon.get())
        lat2, lon2 = float(self.max_lat.get()), float(self.max_lon.get())
        dep1, dep2 = float(self.min_dep.get()), float(self.max_dep.get())
        if self.cat_path.cget("state") != "disabled":
            # acquire station info
            st = self.staz_path.get()
            delim = self.delm.get()
            cod = int(self.code_col.get())
            lat = int(self.lat_col.get())
            lat_f = self.lat_form.get()
            lon = int(self.lon_col.get())
            lon_f = self.lon_form.get()
            dep = int(self.elev_col.get())
            dep_f = self.elev_form.get()
            staz_df = pd.DataFrame(columns=['Name', 'Longitude', 'Latitude', 'Elevation'])
            with open(st, errors='ignore') as o:
                f_st = o.readlines()
            for s in f_st:
                s = s.split(delim)
                staz_df.loc[len(staz_df)] = [s[cod],float(s[lon]),float(s[lat]),float(s[dep])]
            self.staz_df = staz_df.copy()
            convert_depth(dep_f, staz_df, 'Elevation')
            convert_coord(lat_f, staz_df, 'Latitude')
            convert_coord(lon_f, staz_df, 'Longitude')
            # acquire catalog info
            cat = self.cat_path.get()
            form = self.cat_format.get()
            if self.radio_var.get() == 0:
                ev_list = [f for f in glob.glob(f'{cat}/{self.obspy_formats[form]}')]
                if form == "Hypoellipse (.out)":
                    ev_cat = Catalog()
                    for e  in ev_list:
                        ev = parse_hypoellipse(e)
                        for e_v in ev:
                            ev_cat.append(e_v)
                else:
                    ev_cat = Catalog()
                    for e  in ev_list:
                        ev = read_events(e)
                        ev_cat.append(ev[0])
            elif self.radio_var.get() == 1:
                ev_cat = read_events(cat)
            # resize staz
            mask_lat = (staz_df['Latitude'] > lat1) & (staz_df['Latitude'] < lat2)
            staz_df = staz_df[mask_lat]
            mask_lon = (staz_df['Longitude'] > lon1) & (staz_df['Longitude'] < lon2)
            staz_df = staz_df[mask_lon]
            mask_dep = (staz_df['Elevation'] > dep1) & (staz_df['Elevation'] < dep2)
            self.staz_df_r = staz_df[mask_dep]
            # resize earthquake
            ev_copy_cat = Catalog()
            for eq in ev_cat:
                if eq.preferred_origin() == None:
                    origin = eq.origins[0]
                else:
                    origin = eq.preferred_origin()
                or_lon = origin.longitude
                or_lat = origin.latitude
                or_dep = origin.depth/1000
                if lon1 < or_lon < lon2 and lat1 < or_lat < lat2 and \
                    dep1 < or_dep < dep2:
                    ev_copy_cat.append(eq)
            self.ev_cat = ev_copy_cat.copy()

    def get_velocity(self):
        import matplotlib.pyplot as plt
        try:
            max_depth = float(self.max_dep.get())
            min_depth = float(self.min_dep.get())
        except:
            max_depth = 600
            min_depth = 0
        if self.radio_var2.get() == 1:
            self.gradual_velocity = 1
            try:
                self.min_VP = float(self.min_velP.get())
                self.max_VP = float(self.max_velP.get())
                self.min_VS = float(self.min_velS.get())
                self.max_VS = float(self.max_velS.get())
                self.dep_range = max_depth - min_depth
            except:
                # data from IASP91
                self.min_VP = 5.8
                self.max_VP = 9.9984
                self.min_VS = 3.36
                self.max_VS = 5.4728
                self.dep_range = max_depth - min_depth
        elif self.radio_var2.get() == 0:
            self.velocity_model = self.vel_file.get()
            self.gradual_velocity = 0
            delm_list = [";"," ", "|"]
            try:
                with open(self.velocity_model) as f_v:
                    file_v = f_v.readlines()
                while '\n' in file_v:
                    file_v.remove('\n')
                self.vel_mod = pd.DataFrame(columns=['Depth','Pvel','Svel'])
                v_line = file_v[0]
                for delm in delm_list:
                    if len(v_line.split(delm)) > 1:
                        break
                for v_line in file_v:
                    v_line = v_line.split(delm)
                    while '' in v_line:
                        v_line.remove('')
                    while '\n' in v_line:
                        v_line.remove('\n')
                    if v_line != []:
                        if v_line[-1][-1]=='\n':
                            self.vel_mod.loc[len(self.vel_mod)] = [float(v_line[0]), float(v_line[1]),
                                             float(v_line[-1][:-1])]
                        else:
                            self.vel_mod.loc[len(self.vel_mod)] = [float(v_line[0]), float(v_line[1]),
                                             float(v_line[-1])]
                mx_dp = self.vel_mod.loc[np.searchsorted(\
                            np.array(self.vel_mod.Depth),float(self.max_dep.get()))]
                mi_dp = self.vel_mod.loc[np.searchsorted(\
                            np.array(self.vel_mod.Depth),float(self.min_dep.get()))]
                self.min_VP = mi_dp.Pvel
                self.max_VP = mx_dp.Pvel
                self.min_VS = mi_dp.Svel
                self.max_VS = mx_dp.Svel
                self.dep_range = mx_dp.Depth-mi_dp.Depth
            except:
                print("Invalid velocity file, using IASP91")
                self.gradual_velocity = 1
                self.min_VP = 5.8
                self.max_VP = 9.9984
                self.min_VS = 3.36
                self.max_VS = 5.4728
                self.dep_range = max_depth - min_depth
        plot_path = os.getcwd() + '/cat_plots/'
        plot_path = plot_path.replace('\\','/')
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)
        if self.gradual_velocity == 0:
            plt.figure()
            plt.step(self.vel_mod['Pvel'], self.vel_mod['Depth'], 
                     where='post', label='Vp', color='blue', linewidth=2)
            plt.step(self.vel_mod['Svel'], self.vel_mod['Depth'], 
                     where='post', label='Vs', color='red', linewidth=2)
            plt.gca().invert_yaxis()
            plt.xlabel('Velocity (km/s)')
            plt.ylabel('Depth (km)')
            plt.title('1D velocity profile')
            plt.legend()
            plt.grid()
            plt.show()
            plt.savefig(plot_path+'Velocity_profile.png')
            plt.savefig(plot_path+'Velocity_profile.pdf')
            plt.savefig(plot_path+'Velocity_profile.svg')
        elif self.gradual_velocity == 1:
            plt.figure()
            depths = range(int(min_depth), int(max_depth))
            coeff =  (self.max_VP-self.min_VP)/self.dep_range
            plt.plot((self.min_VP + coeff*depths), depths,
                          label='Vp', color='blue', linewidth=2)
            coeff =  (self.max_VS-self.min_VS)/self.dep_range
            plt.plot((self.min_VS + coeff*depths), depths,
                          label='Vs', color='red', linewidth=2)
            plt.gca().invert_yaxis()
            plt.xlabel('Velocity (km/s)')
            plt.ylabel('Depth (km)')
            plt.title('1D velocity profile')
            plt.legend()
            plt.grid()
            plt.show()
            plt.savefig(plot_path+'Velocity_profile.png')
            plt.savefig(plot_path+'Velocity_profile.pdf')
            plt.savefig(plot_path+'Velocity_profile.svg')

    def perform_raytracing(self):
        import ttcrpy.rgrid as rg
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.ticker import FuncFormatter
        import matplotlib.gridspec as gridspec
        from mpl_toolkits.basemap import Basemap
        # get geometry
        self.resize_staz_eqs()
        self.get_velocity()
        y_min, x_min = float(self.min_lat.get()), float(self.min_lon.get())
        y_max, x_max = float(self.max_lat.get()), float(self.max_lon.get())
        z_min, z_max = float(self.min_dep.get()), float(self.max_dep.get())
        if abs(x_max - x_min) > 50 or abs(y_max - y_min) > 50 or abs(z_max - z_min) > 100:
            step = 10
        elif abs(x_max - x_min) < 5 or abs(y_max - y_min) < 5 or abs(z_max - z_min) < 10:
            step = 0.5
        elif abs(x_max - x_min) < 2 or abs(y_max - y_min) < 2 or abs(z_max - z_min) < 5:
            step = 0.1
        else:
            step = 1
        # build grid
        x, y, z = np.arange(start=x_min-step, stop=x_max+step, step=step, dtype=float),\
                    np.arange(start=y_min-step, stop=y_max+step, step=step, dtype=float),\
                    np.arange(start=z_min-step, stop=z_max+step, step=step, dtype=float)
        coeff =  (self.max_VP-self.min_VP)/self.dep_range
        V = np.empty((x.size, y.size, z.size), dtype=float)
        for n in range(z.size):
            V[:, :, n] = self.min_VP + coeff*z[n]
        hitmap = np.zeros(V.shape, dtype=int)
        grid = rg.Grid3d(x, y, z, cell_slowness=False)
        slowness = 1./V
        stz_sel = pd.DataFrame(columns=list(self.staz_df_r.columns))
        ev_sel = pd.DataFrame(columns=["Longitude", "Latitude", "Depth"])
        # calculating rays
        for event in self.ev_cat:
            if event.preferred_origin() == None:
                origin = event.origins[0]
            else:
                origin = event.preferred_origin()
            or_lon = origin.longitude
            or_lat = origin.latitude
            or_dep = origin.depth/1000
            ev_sel.loc[len(ev_sel)] = [or_lon, or_lat, or_dep]
            src = np.array([[or_lon, or_lat, or_dep]])
            rcv = []
            for pick in event.picks:
                if pick.phase_hint == 'P':
                    staz = pick.waveform_id.station_code
                    station = self.staz_df_r[self.staz_df_r["Name"] == staz.upper()]
                    if station.empty:
                        station2 = self.staz_df[self.staz_df["Name"] == staz.upper()]
                        if station2.empty:
                            print(f"{staz.upper()} not found in the list")
                        continue
                    slon, slat, sdep = float(station['Longitude'].iloc[0]),\
                    float(station['Latitude'].iloc[0]),\
                    float(station['Elevation'].iloc[0])
                    stz_sel = pd.concat([stz_sel, station])
                    rcv.append([slon, slat, sdep])
            if rcv == []:
                print("Not station found for the current source")
                continue
            rcv = np.array(rcv)
            try:
                _, rays = grid.raytrace(src, rcv, slowness, return_rays=True)
            except:
                print(f"Error in raytracing for\nssrc: {src}\nrcv: {rcv}")
            mx_ray_dp = 0
            for ray in rays:
                for cell in ray:
                    hitmap[np.argmin(np.abs(x - cell[0])),
                           np.argmin(np.abs(y - cell[1])),
                           np.argmin(np.abs(z - cell[2]))] += 1
                    if float(cell[2]) > mx_ray_dp:
                        mx_ray_dp = float(cell[2])
        # plot hitmaps
        plot_path = os.getcwd() + '/cat_plots/'
        plot_path = plot_path.replace('\\','/')
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)
        stz_sel = stz_sel.drop_duplicates()
        hitmap_normalized = hitmap / np.max(hitmap)
        colors = cm.Greens_r(hitmap_normalized)
        fig = plt.figure()
        gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1])
        ax = fig.add_subplot(gs[0, 0], projection='3d')
        ax.voxels(hitmap > 0, facecolors=colors, alpha=0.3)
        ax.scatter((stz_sel["Longitude"]-(x_min-step))/step,
                   (stz_sel["Latitude"]-(y_min-step))/step,
                   -stz_sel["Elevation"]/1000, marker='^', c='blue')
        ax.scatter((ev_sel["Longitude"]-(x_min-step))/step,
                   (ev_sel["Latitude"]-(y_min-step))/step,
                   ev_sel["Depth"], marker='*', c='red')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_zlabel('Depth')
        ax.set_zlim(z_max, z_min)
        ax.set_xlim(0, x.size)
        ax.set_ylim(0, y.size)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda val, _: f"{val + x_min:.1f}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda val, _: f"{val + y_min:.1f}"))
        ax0 = fig.add_subplot(gs[0, 1])
        h_section = np.sum(hitmap, axis=2)/np.max(hitmap)
        im0 = ax0.imshow(h_section.T, extent=(x.min(), x.max(), y.min(), y.max()),
            origin='lower', cmap='Greens', aspect='auto', vmin=0, vmax=1)
        ax0.set_title("Horizontal section")
        ax0.set_xlabel("Longitude")
        ax0.set_ylabel("Latitude")
        ax0.scatter(stz_sel["Longitude"], stz_sel["Latitude"],
                   marker='^', c='blue', alpha=0.5)
        fig.colorbar(im0, ax=ax0, label=f'Ray norm.\n(max {int(np.max(hitmap))} rays)')
        ax1 = fig.add_subplot(gs[1, 0])
        ns_section = np.sum(hitmap, axis=0)/np.max(hitmap)
        im1 = ax1.imshow(ns_section.T, extent=(y.min(), y.max(), z.min(), z.max()),
            origin='lower', cmap='Greens', aspect='auto', vmin=0, vmax=1)
        ax1.set_title("N-S vertical section")
        ax1.set_xlabel("Latitude")
        ax1.set_ylabel("Depth")
        ax1.invert_yaxis()
        fig.colorbar(im1, ax=ax1, label=f'Ray norm.\n(max {int(np.max(hitmap))} rays)')
        ax2 = fig.add_subplot(gs[1, 1])
        ew_section = np.sum(hitmap, axis=1)/np.max(hitmap)
        im2 = ax2.imshow(ew_section.T, extent=(x.min(), x.max(), z.min(), z.max()),
            origin='lower', cmap='Greens', aspect='auto', vmin=0, vmax=1)
        ax2.set_title("W-E vertical section")
        ax2.set_xlabel("Longitude")
        ax2.set_ylabel("Depth")
        ax2.invert_yaxis()
        fig.colorbar(im2, ax=ax2, label=f'Ray norm.\n(max {int(np.max(hitmap))} rays)')
        plt.tight_layout()
        plt.savefig(plot_path+'Hitmaps.png')
        plt.savefig(plot_path+'Hitmaps.pdf')
        plt.savefig(plot_path+'Hitmaps.svg')
        # plot earthquakes and stations
        m = Basemap(llcrnrlon=x_min - step,
                llcrnrlat=y_min - step,
                urcrnrlon=x_max + step,
                urcrnrlat=y_max + step,
                lat_0=(y_max - y_min)/2,
                lon_0=(x_max - x_min)/2,
                projection='merc',
                resolution = 'h',
                area_thresh=10000.,
                )
        m.drawparallels(np.arange(y_min-step,
                                  y_max+step, step),
                        labels=[1, 0, 0, 0], linewidth=0.2)
        m.drawmeridians(np.arange(x_min-step,
                                  x_max+step, step),
                        labels=[0, 0, 0, 1], linewidth=0.2)
        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()
        m.drawmapboundary(fill_color='#46bcec')
        m.fillcontinents(color = 'white',lake_color='#46bcec')
        rm_lons, rm_lats = m(stz_sel['Longitude'].values,
                             stz_sel['Latitude'].values)
        sm_lons, sm_lats = m(ev_sel['Longitude'].values,
                             ev_sel['Latitude'].values)
        m.scatter(rm_lons, rm_lats, marker = '^', color='b', zorder=5)
        m.scatter(sm_lons, sm_lats, marker = 'o', color='r', zorder=1)
        plt.savefig(plot_path+'GIS_map.png')
        plt.savefig(plot_path+'GIS_map.pdf')
        plt.savefig(plot_path+'GIS_map.svg')

    def compile_input_files(self):
        import shutil
        self.fmtomo_path = self.fmtomo_folder.get()
        if self.eqs_switch.get() == 1:
            if self.checked_eqs != 1:
                self.generate_input()
        grid3dgiP = f'''ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Specify number of layers (= number of interfaces -1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1                     c: Number of layers in model
1                     c: Number of velocity grid types (1 or 2)
0.2                   c: Pinchout distance (km) (>=0.0)
vgrids.in             c: Output velocity grid file
interfaces.in         c: Output interface grid file
145678                c: Seed for random number generation (int)
1.5                   c: Minimum permitted velocity (km/s)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set 3-D grid size and location. Note that all layer
c velocity grids have the same spatial dimension, but can
c have different node densities. Interface grids have the
c same node distribution.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
{self.min_dep.get()}       {self.max_dep.get()}       c: Radial range (top-bottom) of grid (km)
{self.min_lat.get()}   {self.max_lat.get()}         c: Latitudinal range (N-S) of grid (degrees)
{self.min_lon.get()}   {self.max_lon.get()}         c: Longitudinal range (E-W) of grid (degrees)
6371.0                c: Earth radius
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set up propagation grid file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
propgrid.in           c: Name of propagation grid file
{self.n_points_depP.get()}    {self.n_points_latP.get()}    {self.n_points_lonP.get()}        c: Number of points in rad lat, long
{self.refined_grid.get()}     {self.n_cell_refined.get()}              c: Refine factor & no. of local cells
0.05                  c: Cushion factor for prop grid (<<1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c First, set up the velocity grids
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1                     c: Checkerboard polarity (1 or -1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set velocity grid values for layer 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
{self.n2_points_depP.get()}         {self.n2_points_depS.get()}          c: Number of radial grid points (type 1 & 2)
{self.n2_points_latP.get()}        {self.n2_points_latS.get()}           c: Number of grid points in theta (N-S)
{self.n2_points_lonP.get()}        {self.n2_points_lonS.get()}          c: Number of grid points in phi (E-W)
{self.radio_var2.get()}                     c: Use model (0) or constant gradient (1)
P                     c: Use P or S velocity model
{self.vel_file.get()}             c: Velocity model (option 0)
1                     c: Dimension of velocity model (1=1-D,3=3-D)
{self.min_velP.get()}       {self.min_velS.get()}         c: Velocity at origin (km/s) (option 1)
{self.max_velP.get()}       {self.max_velS.get()}         c: Velocity at maximum depth (km/s) (option 1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply random structure to layer 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add random structure (0=no,1=yes)
1.5                   c: Standard deviation of Gaussian noise
1                     c: Add a priori model covariance (0=no,1=yes)?
0.3                   c: Diagonal elements of covariance matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply checkerboard to layer 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
{self.check_switch.get()}                     c: Add checkerboard (0=no,1=yes)
{self.pert_value2.get()}                   c: Maximum perturbation of vertices
{self.check_size.get()}                     c: Checkerboard size (NxNxN)
0                     c: Use spacing (0=no,1=yes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally, apply spikes to layer 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
{self.spike_switch.get()}                     c: Apply spikes (0=no,1=yes)
1                     c: Number of spikes
{self.pert_value.get()}                  c: Amplitude of spike 1                  
{self.spike_dep.get()}  {self.spike_lat.get()}  {self.spike_lon.get()}   c: Coordinates of spike 1 (depth, lat, long)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Now, set up the interface grids
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
12                    c: Number of grid points in theta (N-S)
12                    c: Number of grid points in phi (E-W)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set up interface grid for interface 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Obtain grid from external file (0=no,1=yes)
interface1.z          c: External interface grid file (option 1)
{float(self.min_dep.get())-1}                   c: Height of NW grid point (option 0)
{float(self.min_dep.get())-1}                   c: Height of NE grid point (option 0)
{float(self.min_dep.get())-1}                   c: Height of SW grid point (option 0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply random structure to interface 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add random structure (0=no,1=yes)
2.00                  c: Standard deviation of Gaussian noise
1                     c: Add a priori model covariance (0=no,1=yes)?
0.3                   c: Diagonal elements of covariance matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply checkerboard to interface 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add checkerboard (0=no,1=yes)
5.00                  c: Maximum perturbation of vertices
2                     c: Checkerboard size (NxN)
1                     c: Use spacing (0=no,1=yes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally, apply spikes to interface 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Apply spikes (0=no,1=yes)
2                     c: Number of spikes
-5.00                 c: Amplitude of spike 1
-42.1  146.4          c: Coordinates of spike 1
5.00                  c: Amplitude of spike 2
-41.4  146.4          c: Coordinates of spike 2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set up interface grid for interface 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Obtain grid from external file (0=no,1=yes)
interface2.z          c: External interface grid file (option 1)
{float(self.max_dep.get())+1}                 c: Height of NW grid point (option 0)
{float(self.max_dep.get())+1}                 c: Height of NE grid point (option 0)
{float(self.max_dep.get())+1}                 c: Height of SW grid point (option 0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply random structure to interface 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add random structure (0=no,1=yes)
5.0                   c: Standard deviation of Gaussian noise
1                     c: Add a priori model covariance (0=no,1=yes)?
3.0                   c: Diagonal elements of covariance matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply checkerboard to interface 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add checkerboard (0=no,1=yes)
10.0                  c: Maximum perturbation of vertices
1                     c: Checkerboard size (NxN)
1                     c: Use spacing (0=no,1=yes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally, apply spikes to interface 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Apply spikes (0=no,1=yes)
2                     c: Number of spikes
-9.00                 c: Amplitude of spike 1
-32.2  140.0          c: Coordinates of spike 1
9.00                  c: Amplitude of spike 2
-33.4  142.0          c: Coordinates of spike 2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
'''
        grid3dgiS = f'''ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Specify number of layers (= number of interfaces -1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1                     c: Number of layers in model
1                     c: Number of velocity grid types (1 or 2)
0.2                   c: Pinchout distance (km) (>=0.0)
vgrids.in             c: Output velocity grid file
interfaces.in         c: Output interface grid file
145678                c: Seed for random number generation (int)
1.5                   c: Minimum permitted velocity (km/s)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set 3-D grid size and location. Note that all layer
c velocity grids have the same spatial dimension, but can
c have different node densities. Interface grids have the
c same node distribution.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
{self.min_dep.get()}       {self.max_dep.get()}       c: Radial range (top-bottom) of grid (km)
{self.min_lat.get()}   {self.max_lat.get()}         c: Latitudinal range (N-S) of grid (degrees)
{self.min_lon.get()}   {self.max_lon.get()}         c: Longitudinal range (E-W) of grid (degrees)
6371.0                c: Earth radius
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set up propagation grid file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
propgrid.in           c: Name of propagation grid file
{self.n_points_depS.get()}    {self.n_points_latS.get()}    {self.n_points_lonS.get()}        c: Number of points in rad lat, long
{self.refined_grid.get()}     {self.n_cell_refined.get()}              c: Refine factor & no. of local cells
0.05                  c: Cushion factor for prop grid (<<1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c First, set up the velocity grids
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1                     c: Checkerboard polarity (1 or -1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set velocity grid values for layer 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
{self.n2_points_depP.get()}         {self.n2_points_depS.get()}          c: Number of radial grid points (type 1 & 2)
{self.n2_points_latP.get()}        {self.n2_points_latS.get()}           c: Number of grid points in theta (N-S)
{self.n2_points_lonP.get()}        {self.n2_points_lonS.get()}          c: Number of grid points in phi (E-W)
{self.radio_var2.get()}                     c: Use model (0) or constant gradient (1)
S                     c: Use P or S velocity model
{self.vel_file.get()}             c: Velocity model (option 0)
1                     c: Dimension of velocity model (1=1-D,3=3-D)
{self.min_velP.get()}       {self.min_velS.get()}         c: Velocity at origin (km/s) (option 1)
{self.max_velP.get()}       {self.max_velS.get()}         c: Velocity at maximum depth (km/s) (option 1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply random structure to layer 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add random structure (0=no,1=yes)
1.5                   c: Standard deviation of Gaussian noise
1                     c: Add a priori model covariance (0=no,1=yes)?
0.3                   c: Diagonal elements of covariance matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply checkerboard to layer 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
{self.check_switch.get()}                     c: Add checkerboard (0=no,1=yes)
{self.pert_value2.get()}                   c: Maximum perturbation of vertices
{self.check_size.get()}                     c: Checkerboard size (NxNxN)
0                     c: Use spacing (0=no,1=yes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally, apply spikes to layer 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
{self.spike_switch.get()}                     c: Apply spikes (0=no,1=yes)
1                     c: Number of spikes
{self.pert_value.get()}                  c: Amplitude of spike 1                  
{self.spike_dep.get()}  {self.spike_lat.get()}  {self.spike_lon.get()}   c: Coordinates of spike 1 (depth, lat, long)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Now, set up the interface grids
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
12                    c: Number of grid points in theta (N-S)
12                    c: Number of grid points in phi (E-W)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set up interface grid for interface 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Obtain grid from external file (0=no,1=yes)
interface1.z          c: External interface grid file (option 1)
{float(self.min_dep.get())-1}                   c: Height of NW grid point (option 0)
{float(self.min_dep.get())-1}                   c: Height of NE grid point (option 0)
{float(self.min_dep.get())-1}                   c: Height of SW grid point (option 0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply random structure to interface 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add random structure (0=no,1=yes)
2.00                  c: Standard deviation of Gaussian noise
1                     c: Add a priori model covariance (0=no,1=yes)?
0.3                   c: Diagonal elements of covariance matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply checkerboard to interface 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add checkerboard (0=no,1=yes)
5.00                  c: Maximum perturbation of vertices
2                     c: Checkerboard size (NxN)
1                     c: Use spacing (0=no,1=yes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally, apply spikes to interface 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Apply spikes (0=no,1=yes)
2                     c: Number of spikes
-5.00                 c: Amplitude of spike 1
-42.1  146.4          c: Coordinates of spike 1
5.00                  c: Amplitude of spike 2
-41.4  146.4          c: Coordinates of spike 2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set up interface grid for interface 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Obtain grid from external file (0=no,1=yes)
interface2.z          c: External interface grid file (option 1)
{float(self.max_dep.get())+1}                 c: Height of NW grid point (option 0)
{float(self.max_dep.get())+1}                 c: Height of NE grid point (option 0)
{float(self.max_dep.get())+1}                 c: Height of SW grid point (option 0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply random structure to interface 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add random structure (0=no,1=yes)
5.0                   c: Standard deviation of Gaussian noise
1                     c: Add a priori model covariance (0=no,1=yes)?
3.0                   c: Diagonal elements of covariance matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally apply checkerboard to interface 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Add checkerboard (0=no,1=yes)
10.0                  c: Maximum perturbation of vertices
1                     c: Checkerboard size (NxN)
1                     c: Use spacing (0=no,1=yes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optionally, apply spikes to interface 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
0                     c: Apply spikes (0=no,1=yes)
2                     c: Number of spikes
-9.00                 c: Amplitude of spike 1
-32.2  140.0          c: Coordinates of spike 1
9.00                  c: Amplitude of spike 2
-33.4  142.0          c: Coordinates of spike 2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
'''
        invert3dP = f'''ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This file contains all required input parameters and files for the
c inversion program "invert3d.f90"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
vgridsref.in                c: Reference velocity grid
vgrids.in                   c: Current velocity grid
interfacesref.in            c: Reference interface grid
interfaces.in               c: Current interface grid
sourcesref.in               c: Reference source coords
sources.in                  c: Current source coordinates
stimes.dat                  c: Current source time perturbations
receivers.in                c: Receiver coordinates
otimes.dat                  c: Observed traveltimes
mtimes.dat                  c: Model traveltimes
rtimesnec.dat               c: Reference teleseismic traveltimes
frechet.in                  c: Frechet derivative parameters
frechet.dat                 c: Frechet derivatives
inviter.in                  c: File indicating current inversion step
propgrid.in                 c: File containing propagation grid parameters
0.07                        c: Minimum distance between interfaces (km)
1.5                         c: Minimum permitted velocity (km/s)
0                           c: Remove mean from predicted teleseisms (0=no,1=yes)
1     {self.dampP.get()}    {self.smoothP.get()}            c: Velocity inversion (0=no,1=yes),damp,smooth
0     0.5    0.02           c: Interface inversion (0=no, 1=yes),damp,smooth
0     15      15            c: Source inversion (0=no, 1=yes),damp1,damp2
{self.n_subs.get()}                          c: Subspace dimension (max=50)
1                           c: Global damping factor (epsilon)
0                           c: Apply second derivative smoothing (0=no,1=yes)
0.02                        c: Global smoothing factor (eta)
6371.0                      c: Earth radius in km
'''
        invert3dS = f'''ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This file contains all required input parameters and files for the
c inversion program "invert3d.f90"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
vgridsref.in                c: Reference velocity grid
vgrids.in                   c: Current velocity grid
interfacesref.in            c: Reference interface grid
interfaces.in               c: Current interface grid
sourcesref.in               c: Reference source coords
sources.in                  c: Current source coordinates
stimes.dat                  c: Current source time perturbations
receivers.in                c: Receiver coordinates
otimes.dat                  c: Observed traveltimes
mtimes.dat                  c: Model traveltimes
rtimesnec.dat               c: Reference teleseismic traveltimes
frechet.in                  c: Frechet derivative parameters
frechet.dat                 c: Frechet derivatives
inviter.in                  c: File indicating current inversion step
propgrid.in                 c: File containing propagation grid parameters
0.07                        c: Minimum distance between interfaces (km)
1.5                         c: Minimum permitted velocity (km/s)
0                           c: Remove mean from predicted teleseisms (0=no,1=yes)
1     {self.dampS.get()}    {self.smoothS.get()}            c: Velocity inversion (0=no,1=yes),damp,smooth
0     0.5    0.02           c: Interface inversion (0=no, 1=yes),damp,smooth
0     15      15            c: Source inversion (0=no, 1=yes),damp1,damp2
{self.n_subs.get()}                          c: Subspace dimension (max=50)
1                           c: Global damping factor (epsilon)
0                           c: Apply second derivative smoothing (0=no,1=yes)
0.02                        c: Global smoothing factor (eta)
6371.0                      c: Earth radius in km
'''
        with open(f"{self.fmtomo_path}/invert_p/mkdir/grid3dgi.in", 'w') as o:
            o.write(grid3dgiP)
        with open(f"{self.fmtomo_path}/invert_s/mkdir/grid3dgi.in", 'w') as o:
            o.write(grid3dgiS)
        with open(f"{self.fmtomo_path}/invert_p/invert3d.in", 'w') as o:
            o.write(invert3dP)
        with open(f"{self.fmtomo_path}/invert_s/mkdir/invert3d.in", 'w') as o:
            o.write(invert3dS)

    def reset_GUI(self, entries, buttons, comboboxes):
        from customtkinter import ThemeManager
        for entry in entries:
            if entry.cget("state") == "disabled":
                entry.configure(state="normal")
                entry.delete(0, 'end')
                entry.configure(state="disabled")
            else:
                entry.delete(0, 'end')
        color = ThemeManager.theme["CTkButton"]["fg_color"]
        hv_color = ThemeManager.theme["CTkButton"]["hover_color"]
        for button in buttons:
            if button.cget("state") == "disabled":
                button.configure(state="normal")
                button.configure(fg_color=color, hover_color=hv_color)
                button.configure(state="disabled")
            else:
                button.configure(fg_color=color, hover_color=hv_color)
        for cbx in comboboxes:
            if cbx.cget("state") == "disabled":
                cbx.configure(state="normal")
                cbx[0].set(cbx[1])
                cbx.configure(state="disabled")
            else:
                cbx[0].set(cbx[1])

## Running GUI
if __name__ == "__main__":
    app = GUI()
    app.mainloop()