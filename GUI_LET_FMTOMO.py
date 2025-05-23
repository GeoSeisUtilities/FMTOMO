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
        self.empty_folder = ctk.CTkEntry(self, placeholder_text="Empty files folder", width=width/4)
        self.empty_folder.configure(state="disabled")
        self.empty_folder.place(x=(width/2), y=offset)
        browse_button_2 = ctk.CTkButton(self, text="Browse", width=width/6, state="disabled",
                            command=lambda: self.browse_folder(self.empty_folder))
        browse_button_2.place(x=((width/2) + offset + width/4), y=offset)

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
        self.browse_cat = ctk.CTkButton(Eqs_data_tab, text="Browse", width=width/6,
                            state="disabled", command=lambda: self.browse_folder(self.fmtomo_folder))
        self.browse_cat.place(x=(offset + 3*width/4), y=6*offset)
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
                         "Custom format": ""}
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

        self.gen_input = ctk.CTkButton(Eqs_data_tab, text="Generate input files", width=(width/2-2*offset),
                             state='disabled', command=lambda: self.generate_input())
        self.gen_input.place(x=offset, y=25*offset)
        eqs_widgets = [self.cat_path, self.is_cat_folder, self.is_cat_file, self.cat_format, self.staz_path,
                       self.browse_cat, self.browse_staz, self.delm, self.code_col, self.lat_col, self.lat_form,
                       self.lon_col, self.lon_form, self.elev_col, self.elev_form, self.PS, self.gen_input]
        custom_format_widgets = [self.delm_cat, self.time_cat_col, self.time_cat_form,
                    self.lat_cat_col, self.lat_cat_form, self.lon_cat_col, self.lon_cat_form,
                    self.elev_cat_col, self.elev_cat_form, self.delm_pick, self.pick_form]

        # direct problem tab
        dir_switch = ctk.CTkSwitch(master=Dir_prob_tab, switch_height=10,
                            text="I need to define the model",
                            command=lambda: self.switch_activation([self.empty_folder, browse_button_2]))
        dir_switch.place(x=10, y=10)


        # inverse problem tab
        inv_switch = ctk.CTkSwitch(master=Inv_prob_tab, switch_height=10,
                            text="I need to set the inversion",
                            command=lambda: self.switch_activation([self.empty_folder, browse_button_2]))
        inv_switch.place(x=10, y=10)

    def browse_folder(self, entry):
        from tkinter import filedialog
        filename = filedialog.askdirectory(initialdir = "/home/", title = "Browse folders",)
        entry.delete(0,)
        entry.insert(0,str(filename))

    def browse_file(self, entry):
        from tkinter import filedialog
        filename=filedialog.askopenfilename(initialdir = "/home/",title = "Open file",
                                            filetypes = (("Text files","*.txt"),
                                            ("dat files","*.dat"),("csv files","*.csv"),
                                            ("All files","*.*")))
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

    def second_switch_activator (self, entry, value, widgets):
        if entry == value:
            self.switch_activation(widgets)
        else:
            state = widgets[0].cget("state")
            if state != "disabled":
                self.switch_activation(widgets)

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
            if self.radio_var == 0:
                ev_list = [f for f in glob.glob(f'{cat}/{self.obspy_formats[form]}')]
                ev_cat = Catalog()
                for e  in ev_list:
                    ev = read_events(e)
                    ev_cat.append(ev[0])
            elif self.radio_var == 1:
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

## Running GUI
if __name__ == "__main__":
    app = GUI()
    app.mainloop()