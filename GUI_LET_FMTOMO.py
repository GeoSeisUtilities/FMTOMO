#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 11:20:28 2025

@author: GeoSeisUtilities
"""

### This is the GUI for performing LET with FMTOMO-Extras

## Importing module
import glob
from obspy.core.event import read_events
from obspy import Catalog
import pandas as pd
import pathlib
import numpy as np
import customtkinter as ctk
from tkinter import IntVar

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
        self.fmtomo_folder = ctk.CTkEntry(self, placeholder_text="FMTOMO folder", width=width/4)
        self.fmtomo_folder.place(x=offset, y=offset)
        browse_button_1 = ctk.CTkButton(self, text="Browse", width=width/6,
                            command=lambda: self.browse_folder(self.fmtomo_folder))
        browse_button_1.place(x=(2*offset + width/4), y=offset)

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
                            state="disabled", command=lambda: self.browse_folder(self.fmtomo_folder)\
                            if self.radio_var.get() == 0 else self.browse_file(self.fmtomo_folder,\
                            list(self.obspy_formats.items())))
        self.browse_cat.place(x=(offset + 3*width/4), y=6*offset)
        self.cat_format = ctk.CTkComboBox(Eqs_data_tab, width=width/4,
                                     values=list(self.obspy_formats.keys()),
                                     command=lambda x: self.second_switch_activator(\
                                     x, 'Custom format', custom_format_widgets))
        self.cat_format.set("Catalog format")
        self.cat_format.place(x=width/2, y=9*offset)
        self.cat_format.configure(state="disabled")
        self.PS_choiche = ['P', 'S', 'P and S']
        self.PS = ctk.CTkComboBox(Eqs_data_tab, width=width/6,
                                     values=self.PS_choiche)
        self.PS.set("P or S")
        self.PS.place(x=(3*width/4 + offset), y=9*offset)
        self.PS.configure(state="disabled")
        self.staz_path = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Station file",
                                width=width/4)
        self.staz_path.place(x=offset, y=3*offset)
        self.staz_path.configure(state="disabled")
        self.browse_staz = ctk.CTkButton(Eqs_data_tab, text="Browse", width=width/6,
                            state="disabled",  command=lambda: self.browse_file(self.staz_path))
        self.browse_staz.place(x=(2*offset + width/4), y=3*offset)
        # station entries
        self.delm = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Delimiter",
                                width=width/4)
        self.delm.place(x=offset, y=6*offset)
        self.delm.configure(state="disabled")
        self.code_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Name (code) column",
                                width=width/4)
        self.code_col.place(x=offset, y=9*offset)
        self.code_col.configure(state="disabled")
        self.coord_format = ["dd mm ss", "dd mm.mm", "ddLmm.mm", "dd-mm-ss",
                        "dd L mm ss", "ddLmmss.sss", "dd.dd"]
        self.lat_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Latitude column",
                                width=width/4)
        self.lat_col.place(x=offset, y=12*offset)
        self.lat_col.configure(state="disabled")
        self.lat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/6,
                                     values=self.coord_format)
        self.lat_form.set("Format")
        self.lat_form.place(x=(width/4 + 2*offset), y=12*offset)
        self.lat_form.configure(state="disabled")
        self.lon_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Longitude column",
                                width=width/4)
        self.lon_col.place(x=offset, y=15*offset)
        self.lon_col.configure(state="disabled")
        self.lon_form = ctk.CTkComboBox(Eqs_data_tab, width=width/6,
                                     values=self.coord_format)
        self.lon_form.set("Format")
        self.lon_form.place(x=(width/4 + 2*offset), y=15*offset)
        self.lon_form.configure(state="disabled")
        self.depth_format = ["Elevation (m)", "Depth (m)", "Elevation (Km)", "Depth (Km)"]
        self.elev_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Elevation column",
                                width=width/4)
        self.elev_col.place(x=offset, y=18*offset)
        self.elev_col.configure(state="disabled")
        self.elev_form = ctk.CTkComboBox(Eqs_data_tab, width=width/6,
                                     values=self.depth_format)
        self.elev_form.set("Format")
        self.elev_form.place(x=(width/4 + 2*offset), y=18*offset)
        self.elev_form.configure(state="disabled")
        # catalog entries
        self.delm_cat = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Cat. Delim.",
                                width=(width/8 - offset))
        self.delm_cat.place(x=width/2, y=13*offset)
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
        self.time_cat_col.configure(state="disabled")
        self.time_cat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/8,
                                     values=list(date_format.keys()))
        self.time_cat_form.set("Format")
        self.time_cat_form.place(x=(5*width/8), y=16*offset)
        self.time_cat_form.configure(state="disabled")
        self.lat_cat_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Lat. col.",
                                width=(width/8 - offset))
        self.lat_cat_col.place(x=width/2, y=19*offset)
        self.lat_cat_col.configure(state="disabled")
        self.lat_cat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/8,
                                     values=self.coord_format)
        self.lat_cat_form.set("Format")
        self.lat_cat_form.place(x=(5*width/8), y=19*offset)
        self.lat_cat_form.configure(state="disabled")
        self.lon_cat_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Lon. col.",
                                width=(width/8 - offset))
        self.lon_cat_col.place(x=width/2, y=22*offset)
        self.lon_cat_col.configure(state="disabled")
        self.lon_cat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/8,
                                     values=self.coord_format)
        self.lon_cat_form.set("Format")
        self.lon_cat_form.place(x=(5*width/8), y=22*offset)
        self.lon_cat_form.configure(state="disabled")
        self.elev_cat_col = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Dep. col.",
                                width=(width/8 - offset))
        self.elev_cat_col.place(x=width/2, y=25*offset)
        self.elev_cat_col.configure(state="disabled")
        self.elev_cat_form = ctk.CTkComboBox(Eqs_data_tab, width=width/8,
                                     values=self.depth_format)
        self.elev_cat_form.set("Format")
        self.elev_cat_form.place(x=(5*width/8), y=25*offset)
        self.elev_cat_form.configure(state="disabled")
        # picks entries
        self.delm_pick = ctk.CTkEntry(Eqs_data_tab, placeholder_text="Pick Delim.",
                                width=(width/8 - offset))
        self.delm_pick.place(x=(3*width/4 + offset), y=13*offset)
        self.delm_pick.configure(state="disabled")
        picking_format = ["Still working on...", "Still working on..."]
        self.pick_form = ctk.CTkComboBox(Eqs_data_tab, width=width/9,
                                     values=picking_format)
        self.pick_form.set("Format")
        self.pick_form.place(x=(3*width/4 + offset), y=16*offset)
        self.pick_form.configure(state="disabled")
        # button
        self.gen_input = ctk.CTkButton(Eqs_data_tab, text="Generate input files", width=(2*(width/5) + 2*offset),
                             state='disabled', command=lambda: self.generate_input())
        self.gen_input.place(x=offset, y=25*offset)
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
        self.min_lat.configure(state="disabled")
        self.max_lat = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Maximum Latitude",
                                width=(width/5 - 2*offset))
        self.max_lat.place(x=(width/5), y=3*offset)
        self.max_lat.configure(state="disabled")
        P_lab = ctk.CTkLabel(Dir_prob_tab, text="Nodes P model")
        P_lab.place(x=(2*(width/5) - offset), y=0)
        S_lab = ctk.CTkLabel(Dir_prob_tab, text="Nodes S model")
        S_lab.place(x=(3*(width/5) - offset), y=0)
        self.n_points_latP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_latP.place(x=(2*(width/5) - offset), y=3*offset)
        self.n_points_latP.configure(state="disabled")
        self.n2_points_latP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_latP.place(x=(2.5*(width/5) - offset), y=3*offset)
        self.n2_points_latP.configure(state="disabled")
        self.n_points_latS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_latS.place(x=(3*(width/5) - offset), y=3*offset)
        self.n_points_latS.configure(state="disabled")
        self.n2_points_latS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_latS.place(x=(3.5*(width/5) - offset), y=3*offset)
        self.n2_points_latS.configure(state="disabled")
        self.min_lon = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Min. Longitude",
                                width=(width/5 - 2*offset))
        self.min_lon.place(x=offset, y=6*offset)
        self.min_lon.configure(state="disabled")
        self.max_lon = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Max. Longitude",
                                width=(width/5 - 2*offset))
        self.max_lon.place(x=(width/5), y=6*offset)
        self.max_lon.configure(state="disabled")
        self.n_points_lonP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_lonP.place(x=(2*(width/5) - offset), y=6*offset)
        self.n_points_lonP.configure(state="disabled")
        self.n2_points_lonP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_lonP.place(x=(2.5*(width/5) - offset), y=6*offset)
        self.n2_points_lonP.configure(state="disabled")
        self.n_points_lonS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_lonS.place(x=(3*(width/5) - offset), y=6*offset)
        self.n_points_lonS.configure(state="disabled")
        self.n2_points_lonS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_lonS.place(x=(3.5*(width/5) - offset), y=6*offset)
        self.n2_points_lonS.configure(state="disabled")
        self.min_dep = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Min. Depth (>0)",
                                width=(width/5 - 2*offset))
        self.min_dep.place(x=offset, y=9*offset)
        self.min_dep.configure(state="disabled")
        self.max_dep = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Max. Depth (<0)",
                                width=(width/5 - 2*offset))
        self.max_dep.place(x=(width/5), y=9*offset)
        self.max_dep.configure(state="disabled")
        self.n_points_depP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_depP.place(x=(2*(width/5) - offset), y=9*offset)
        self.n_points_depP.configure(state="disabled")
        self.n2_points_depP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_depP.place(x=(2.5*(width/5) - offset), y=9*offset)
        self.n2_points_depP.configure(state="disabled")
        self.n_points_depS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Prop.",
                                width=(width/10 - offset))
        self.n_points_depS.place(x=(3*(width/5) - offset), y=9*offset)
        self.n_points_depS.configure(state="disabled")
        self.n2_points_depS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Vel.",
                                width=(width/10 - offset))
        self.n2_points_depS.place(x=(3.5*(width/5) - offset), y=9*offset)
        self.n2_points_depS.configure(state="disabled")
        self.refined_grid = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Refine grid dim.",
                                width=(width/5 - 3*offset))
        self.refined_grid.configure(state="disabled")
        self.refined_grid.place(x=(4*(width/5) - offset), y=4.5*offset)
        self.n_cell_refined = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Refine grid ext.",
                                width=(width/5 - 3*offset))
        self.n_cell_refined.configure(state="disabled")
        self.n_cell_refined.place(x=(4*(width/5) - offset), y=7.5*offset)
        # estimate spacing
        self.est_dist = ctk.CTkButton(Dir_prob_tab, text="Estimate point spacing", width=(width - 5*offset),
                             state='disabled', command=lambda: self.estimate_spacing())
        self.est_dist.place(x=offset, y=13*offset)
        estimated_prop = ctk.CTkLabel(Dir_prob_tab, text="Propagation grid spacing (km)")
        estimated_prop.place(x=offset, y=16*offset)
        self.est_latP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lat. spac.",
                                width=(width/6 - 2*offset))
        self.est_latP.place(x=offset, y=19*offset)
        self.est_latP.configure(state="disabled")
        self.est_lonP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lon. spac.",
                                width=(width/6 - 2*offset))
        self.est_lonP.place(x=offset, y=22*offset)
        self.est_lonP.configure(state="disabled")
        self.est_depP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Dep. spac.",
                                width=(width/6 - 2*offset))
        self.est_depP.place(x=offset, y=25*offset)
        self.est_depP.configure(state="disabled")
        self.est_latS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lat. spac.",
                                width=(width/6 - 2*offset))
        self.est_latS.place(x=width/6, y=19*offset)
        self.est_latS.configure(state="disabled")
        self.est_lonS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lon. spac.",
                                width=(width/6 - 2*offset))
        self.est_lonS.place(x=width/6, y=22*offset)
        self.est_lonS.configure(state="disabled")
        self.est_depS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Dep. spac.",
                                width=(width/6 - 2*offset))
        self.est_depS.place(x=width/6, y=25*offset)
        self.est_depS.configure(state="disabled")
        estimated_vel = ctk.CTkLabel(Dir_prob_tab, text="Velocity grid spacing (km)")
        estimated_vel.place(x=width/3, y=16*offset)
        self.est_lat2P = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lat. spac.",
                                width=(width/6 - 2*offset))
        self.est_lat2P.place(x=width/3, y=19*offset)
        self.est_lat2P.configure(state="disabled")
        self.est_lon2P = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lon. spac.",
                                width=(width/6 - 2*offset))
        self.est_lon2P.place(x=width/3, y=22*offset)
        self.est_lon2P.configure(state="disabled")
        self.est_dep2P = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Dep. spac.",
                                width=(width/6 - 2*offset))
        self.est_dep2P.place(x=width/3, y=25*offset)
        self.est_dep2P.configure(state="disabled")
        self.est_lat2S = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lat. spac.",
                                width=(width/6 - 2*offset))
        self.est_lat2S.place(x=(width/2 - offset), y=19*offset)
        self.est_lat2S.configure(state="disabled")
        self.est_lon2S = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Lon. spac.",
                                width=(width/6 - 2*offset))
        self.est_lon2S.place(x=(width/2 - offset), y=22*offset)
        self.est_lon2S.configure(state="disabled")
        self.est_dep2S = ctk.CTkEntry(Dir_prob_tab, placeholder_text="Dep. spac.",
                                width=(width/6 - 2*offset))
        self.est_dep2S.place(x=(width/2 - offset), y=25*offset)
        self.est_dep2S.configure(state="disabled")
        tot_n_points = ctk.CTkLabel(Dir_prob_tab, text="Number of propagation nodes")
        tot_n_points.place(x=(2*(width/3) - offset), y=16*offset)
        self.est_tot_pointsP = ctk.CTkEntry(Dir_prob_tab, placeholder_text="P nodes",
                                width=(width/6 - 2*offset))
        self.est_tot_pointsP.configure(state="disabled")
        self.est_tot_pointsP.place(x=(2*(width/3) - offset), y=19*offset)
        self.est_tot_pointsS = ctk.CTkEntry(Dir_prob_tab, placeholder_text="S nodes",
                                width=(width/6 - 2*offset))
        self.est_tot_pointsS.configure(state="disabled")
        self.est_tot_pointsS.place(x=(5*(width/6) - 2*offset), y=19*offset)
        tot_n_points2 = ctk.CTkLabel(Dir_prob_tab, text="Number of velocity nodes")
        tot_n_points2.place(x=2*(width/3) - offset, y=22*offset)
        self.est_tot_points2P = ctk.CTkEntry(Dir_prob_tab, placeholder_text="P nodes",
                                width=(width/6 - 2*offset))
        self.est_tot_points2P.configure(state="disabled")
        self.est_tot_points2P.place(x=(2*(width/3) - offset), y=25*offset)
        self.est_tot_points2S = ctk.CTkEntry(Dir_prob_tab, placeholder_text="S nodes",
                                width=(width/6 - 2*offset))
        self.est_tot_points2S.configure(state="disabled")
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
        self.n_iter.place(x=offset, y=3*offset)
        self.n_subs = ctk.CTkEntry(Inv_prob_tab, placeholder_text="N. of subspaces",
                                width=(width/5))
        self.n_subs.configure(state="disabled")
        self.n_subs.place(x=(width/5 + 2*offset), y=3*offset)
        self.dampP = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Damping value P",
                                width=(width/5))
        self.dampP.configure(state="disabled")
        self.dampP.place(x=offset, y=6*offset)
        self.smoothP = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Smoothing value P",
                                width=(width/5))
        self.smoothP.configure(state="disabled")
        self.smoothP.place(x=offset, y=9*offset)
        self.dampS = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Damping value S",
                                width=(width/5))
        self.dampS.configure(state="disabled")
        self.dampS.place(x=(width/5 + 2*offset), y=6*offset)
        self.smoothS = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Smoothing value S",
                                width=(width/5))
        self.smoothS.configure(state="disabled")
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
        self.vel_file.configure(state="disabled")
        self.vel_filetypes = (("Velocity files","*.vel"),("Text files","*.txt"),
                              ("Dat files","*.dat"),("csv files","*.csv"),("All files","*.*"))
        self.browse_vel = ctk.CTkButton(Inv_prob_tab, text="Browse", width=width/6,
                            state="disabled",  command=lambda: self.browse_file(self.vel_file,
                            self.vel_filetypes))
        self.browse_vel.place(x=(2*offset + 3*(width/4)), y=6*offset)
        self.min_velP = ctk.CTkEntry(Inv_prob_tab, placeholder_text="P velocity top",
                                width=(width/5+offset))
        self.min_velP.place(x=width/2, y=10*offset)
        self.min_velP.configure(state="disabled")
        self.min_velS = ctk.CTkEntry(Inv_prob_tab, placeholder_text="S velocity top",
                                width=(width/5+offset))
        self.min_velS.place(x=(2*offset + 7*(width/10)), y=10*offset)
        self.min_velS.configure(state="disabled")
        self.max_velP = ctk.CTkEntry(Inv_prob_tab, placeholder_text="P velocity bottom",
                                width=(width/5+offset))
        self.max_velP.place(x=width/2, y=13*offset)
        self.max_velP.configure(state="disabled")
        self.max_velS = ctk.CTkEntry(Inv_prob_tab, placeholder_text="S velocity bottom",
                                width=(width/5+offset))
        self.max_velS.place(x=(2*offset + 7*(width/10)), y=13*offset)
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
        self.pert_value.configure(state="disabled")
        self.spike_lat = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Lat",
                                width=(width/6 - 3*offset))
        self.spike_lat.place(x=offset, y=24*offset)
        self.spike_lat.configure(state="disabled")
        self.spike_lon = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Lon",
                                width=(width/6 - 3*offset))
        self.spike_lon.place(x=(width/6 - offset), y=24*offset)
        self.spike_lon.configure(state="disabled")
        self.spike_dep = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Dep",
                                width=(width/6 - 3*offset))
        self.spike_dep.place(x=(width/3 - 3*offset), y=24*offset)
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
        self.pert_value2.configure(state="disabled")
        self.check_size = ctk.CTkEntry(Inv_prob_tab, placeholder_text="Check. dimension",
                                width=(width/5+offset))
        self.check_size.place(x=width/2, y=24*offset)
        self.check_size.configure(state="disabled")
        self.comp_inp = ctk.CTkButton(Inv_prob_tab, text="Compile FMTOMO\n files",
                            width=(width/5 + offset), height=(height/6 - offset),
                            state='disabled', command=lambda: self.compile_input_files())
        self.comp_inp.place(x=(2*offset + 7*(width/10)), y=21*offset)
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
            state = widgets[0].cget("state")
            if state != "disabled":
                self.switch_activation(widgets)

    def third_switch_activator(self, activate, deactivate):
        for wid in activate:
            wid.configure(state="normal")
        for wid2 in deactivate:
            wid2.configure(state="disabled")

    def generate_input(self):
        self.checked_eqs = 1
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
            df = pd.read_csv(st, sep=delim, usecols=[f'col{cod}',f'col{lon}',f'col{lat}',f'col{dep}'], low_memory=True)
            staz_df = df.set_axis(['Name', 'Longitude', 'Latitude', 'Elevation'], axis=1)
            convert_depth(dep_f, staz_df, 'Elevation')
            staz_df['Elevation'] = staz_df['Elevation'].apply(lambda x: -x)
            convert_coord(lat_f, staz_df, 'Latitude')
            convert_coord(lon_f, staz_df, 'Longitude')
            # acquire catalog info
            cat = self.cat_path.get()
            form = self.cat_format.get()
            if self.radio_var.get() == 0:
                ev_list = [f for f in glob.glob(f'{cat}/{self.obspy_formats[form]}')]
                ev_cat = Catalog()
                for e  in ev_list:
                    ev = read_events(e)
                    ev_cat.append(ev[0])
            elif self.radio_var.get() == 1:
                ev_cat = read_events(cat)
            fmtomo_path = self.fmtomo_folder.get()
            ps = self.PS.get()
            if ps == 'P':
                obspy2fmtomo(ev_cat, staz_df, f'{fmtomo_path}/invert_p/mkdir/', ["P"])
            elif ps == 'S':
                obspy2fmtomo(ev_cat, staz_df, f'{fmtomo_path}/invert_s/mkdir/', ["S"])
            elif ps == 'P and S':
                obspy2fmtomo(ev_cat, staz_df, f'{fmtomo_path}/invert_p/mkdir/', ["P"])
                obspy2fmtomo(ev_cat, staz_df, f'{fmtomo_path}/invert_s/mkdir/', ["S"])

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

    def compile_input_files(self):
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

## Running GUI
if __name__ == "__main__":
    app = GUI()
    app.mainloop()