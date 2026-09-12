#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PBUTs SIMULATOR - GUI (Tkinter)
=================================
Interfaccia grafica per il modello Forward Euler (HDF post-diluizione)
delle tossine uremiche legate alle proteine (PBUTs).

Basato su / riutilizza la fisica del modello contenuto in
PBUTs_Excel_Simulator.py (classe PbutsPredictor).

L'utente puo' impostare da interfaccia:
    - Tempo totale di dialisi (min)
    - Tossina analizzata (IS = indoxyl sulfate, pCS = p-cresil solfato)
    - Portata sangue/plasmatica Qp (mL/min)
    - Portata dialisato Qd (mL/min)
    - Filtro utilizzato (Fx100, Solacea21H)
    - Concentrazione plasmatica iniziale di tossina LIBERA (mg/L)
    - Concentrazione plasmatica iniziale di tossina TOTALE (mg/L) [libera + legata]
    - Intervallo di campionamento del grafico (min) - default 5 min
    - (parametri avanzati opzionali: Quf, Qr, dt)

Durante la simulazione il grafico viene aggiornato in tempo reale ogni
intervallo di campionamento impostato (default 5 minuti simulati), cosi'
da vedere l'avanzamento della curva di tossina libera e di tossina totale
mentre la simulazione procede.

Autore: adattato da PBUTs_Excel_Simulator.py
"""

import os
import queue
import threading
import tkinter as tk
from datetime import datetime
from tkinter import ttk, messagebox, filedialog

import numpy as np

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


# =====================================================================
# 1. PARAMETRI DI DEFAULT DEL MODELLO (identici allo script originale)
# =====================================================================

Filter_models = {
    "Fx100": {
        "KoA": 4.08305e-4,   # m^3/min
        "Lfiber": 22.8e-2,   # m
        "N": 16896,          # fibers
        "A": 2.21e-8,        # m^2 fiber area (lato sangue)
        "Ad": 4.30e-8,       # m^2 fiber area (lato dialisato)
        "S_T": 1
    },
    "Solacea21H": {
        "KoA": 5.42956e-4,   # m^3/min
        "Lfiber": 25.4e-2,   # m
        "N": 12960,          # fibers
        "A": 2.78e-8,        # m^2 fiber area (lato sangue)
        "Ad": 4.45e-8,       # m^2 fiber area (lato dialisato)
        "S_T": 1
    }
}

Toxin_parameters = {
    "IS": {   # indoxyl sulfate
        "MwT": 213.21,
        "k1": 3.6e4,
        "G": 0.0285 / 213.21 * 1e-3,
        "Kip_T": 1210e-6,
        "Kic_T": 103e-6,
    },
    "pCS": {  # p-cresol sulfate
        "MwT": 188.22,
        "k1": 1e4,
        "G": 0.0277 / 188.22 * 1e-3,
        "Kip_T": 1060e-6,
        "Kic_T": 98e-6
    }
}

Model_Parameters_default = {
    "n_segments_dialyzer": 40,
    "n_segments_tube": 100,
    "dt": 0.00001,     # min
    "total_time": 240,  # min
    "out_time": 5,      # min - intervallo di campionamento/grafico
    "MwP": 66.1 * 1e3
}

Dialyzation_settings_default = {
    "Qp": 294e-6,        # m^3/min
    "Quf": 95.25862e-6,  # m^3/min
    "Qd": 417.67e-6,     # m^3/min
    "Qr": 82.3276e-6     # m^3/min
}

Patient_default = {
    "Vpl0": 3.775e-3,
    "Vis0": 11.325e-3,
    "Vic": 18.4e-3,
    "Tpl0": 2.357 / 213.21,
    "PTpl0": (31.02 - 2.357) / (213.21),
    "Ppl0": 0.6,
    "Tis0": 2.357 / 213.21,
    "PTis0": (31.02 - 2.357) / (2 * 213.21),
    "Pis0": 0.3,
    "Tic0": 2.357 / 213.21
}

Tube_Parameters_default = {
    "Ltube": 0.5,
    "Atube": 1e-3 ** 2 * 3.14,
}


# =====================================================================
# 2. CLASSE DEL MODELLO (Forward Euler) - identica alla logica originale,
#    con l'aggiunta di un meccanismo di progress-callback per la GUI e
#    l'accumulo in memoria dei risultati per il plotting.
# =====================================================================

class PbutsPredictor:
    def __init__(self, filter_model="Fx100", Toxin="IS", Patient=None, Model_Parameters=None,
                 Tube_Parameters=None, Dialyzation_settings=None, out_name="",
                 out_path=".", progress_callback=None, save_to_file=True):

        Patient = Patient or Patient_default
        Model_Parameters = Model_Parameters or Model_Parameters_default
        Tube_Parameters = Tube_Parameters or Tube_Parameters_default
        Dialyzation_settings = Dialyzation_settings or Dialyzation_settings_default

        for key, value in Filter_models[filter_model].items():
            setattr(self, key, value)
        for key, value in Toxin_parameters[Toxin].items():
            setattr(self, key, value)
        for key, value in Patient.items():
            setattr(self, key, value)
        for key, value in Model_Parameters.items():
            setattr(self, key, value)
        for key, value in Dialyzation_settings.items():
            setattr(self, key, value)
        for key, value in Tube_Parameters.items():
            setattr(self, key, value)

        self.filter_model = filter_model
        self.toxin_name = Toxin
        self.out_path = out_path
        self.out_name = out_name
        self.progress_callback = progress_callback
        self.save_to_file = save_to_file

        self.MwPT = self.MwP + self.MwT
        self.dx = self.Lfiber / self.n_segments_dialyzer
        self.total_time_steps = int(self.total_time / self.dt)
        self.out_steps = max(1, int(self.out_time / self.dt))
        self.dx_tube = self.Ltube / self.n_segments_tube
        self.Qtube = self.Qp - self.Quf + self.Qr

        self.stability_msgs = []
        if self.dt / (self.dx ** 2) > 1:
            self.stability_msgs.append(
                f"Attenzione: condizione di stabilita' NON soddisfatta nel dializzatore "
                f"(dt/dx^2={self.dt / (self.dx ** 2):.2f} > 1). La simulazione potrebbe essere instabile.")
        if self.dt / (self.dx_tube ** 2) > 1:
            self.stability_msgs.append(
                f"Attenzione: condizione di stabilita' NON soddisfatta nel tubo venoso "
                f"(dt/dx^2={self.dt / (self.dx_tube ** 2):.2f} > 1). La simulazione potrebbe essere instabile.")

        self.Pe = self.Quf / self.KoA
        self.Qdx = lambda x: self.Qd + (self.Lfiber - x) / self.Lfiber * self.Quf
        self.Qpx = lambda x: self.Qp - x / self.Lfiber * self.Quf
        self.alpha_t = lambda Vis, Vpl: Vis / (Vis + Vpl)

        self.running = False

        # Risultati campionati (per plotting)
        self.res_time = []
        self.res_T_pl = []
        self.res_PT_pl = []
        self.res_P_pl = []
        self.res_T_is = []
        self.res_T_d_out = []

    def stop(self):
        self.running = False

    def run(self):
        self.running = True
        total_time_steps = self.total_time_steps

        Vpl = np.zeros(total_time_steps)
        Vis = np.zeros(total_time_steps)
        T_pl = np.zeros(total_time_steps)
        PT_pl = np.zeros(total_time_steps)
        P_pl = np.zeros(total_time_steps)
        T_is = np.zeros(total_time_steps)
        PT_is = np.zeros(total_time_steps)
        P_is = np.zeros(total_time_steps)
        T_ic = np.zeros(total_time_steps)
        T_dpl = np.zeros((2, self.n_segments_dialyzer))
        PT_dpl = np.zeros((2, self.n_segments_dialyzer))
        P_dpl = np.zeros((2, self.n_segments_dialyzer))
        T_d = np.zeros((2, self.n_segments_dialyzer))

        T_tpl = np.zeros((2, self.n_segments_tube))
        PT_tpl = np.zeros((2, self.n_segments_tube))
        P_tpl = np.zeros((2, self.n_segments_tube))

        Vpl[0] = self.Vpl0
        Vis[0] = self.Vis0
        Vic = self.Vic
        T_pl[0] = self.Tpl0
        PT_pl[0] = self.PTpl0
        P_pl[0] = self.Ppl0
        T_is[0] = self.Tis0
        PT_is[0] = self.PTis0
        P_is[0] = self.Pis0
        T_ic[0] = self.Tic0
        if self.PTpl0 <= 0:
            if self.progress_callback:
                self.progress_callback(
                    "error", 0, total_time_steps,
                    "La concentrazione totale iniziale deve essere maggiore della concentrazione "
                    "libera iniziale (la quota legata alle proteine non puo' essere nulla o negativa).")
            self.running = False
            return
        self.k2 = self.k1 * (self.Ppl0 * self.Tpl0) / self.PTpl0

        j_dialyzer = np.arange(1, self.n_segments_dialyzer)
        k_dialyzer = j_dialyzer - 1
        j_tube = np.arange(1, self.n_segments_tube)

        start_time = datetime.today().strftime('%Y%m%d_%H%M')
        files = None
        if self.save_to_file:
            os.makedirs(self.out_path, exist_ok=True)
            self.out_results_patient = os.path.join(
                self.out_path, f'FE_results_patient_{self.out_name}_{self.toxin_name}_{start_time}.txt')
            self.out_results_dialyzer = os.path.join(
                self.out_path, f'FE_results_dialyzer_{self.out_name}_{self.toxin_name}_{start_time}.txt')
            patient_file = open(self.out_results_patient, 'w')
            dialyzer_file = open(self.out_results_dialyzer, 'w')
            patient_file.write('Time (min), T_pl (mg/L), PT_pl (mg/L), P_pl (mg/L), '
                                'T_is (mg/L), PT_is (mg/L), P_is (mg/L), T_ic (mg/L)\n')
            dialyzer_file.write('Time (min), T_dpl_out (mg/L), PT_dpl_out (mg/L), '
                                 'P_dpl_out (mg/L), T_d_out (mg/L)\n')
            files = (patient_file, dialyzer_file)

        try:
            for i in range(0, total_time_steps - 1):
                if np.isnan(T_pl[i]):
                    if self.progress_callback:
                        self.progress_callback("error", i, total_time_steps,
                                                "NaN rilevato: la simulazione e' instabile. "
                                                "Prova a ridurre il passo temporale dt.")
                    break
                if not self.running:
                    if self.progress_callback:
                        self.progress_callback("stopped", i, total_time_steps, "Simulazione interrotta dall'utente.")
                    break

                if i % self.out_steps == 0:
                    time_min = i * self.dt
                    self.res_time.append(time_min)
                    self.res_T_pl.append(T_pl[i] * self.MwT)
                    self.res_PT_pl.append(PT_pl[i] * self.MwT)
                    self.res_P_pl.append(P_pl[i] * self.MwP)
                    self.res_T_is.append(T_is[i] * self.MwT)
                    self.res_T_d_out.append(T_d[0, 0] * self.MwT)

                    if files:
                        patient_file, dialyzer_file = files
                        patient_file.write(
                            f'{time_min:.2f}, {T_pl[i] * self.MwT:.5e}, {PT_pl[i] * self.MwT:.5e}, '
                            f'{P_pl[i] * self.MwP:.5e}, {T_is[i] * self.MwT:.5e}, {PT_is[i] * self.MwT:.5e}, '
                            f'{P_is[i] * self.MwP:.5e}, {T_ic[i] * self.MwT:.5e}\n')
                        dialyzer_file.write(
                            f'{time_min:.2f}, {T_dpl[0, -1] * self.MwT:.5e}, {PT_dpl[0, -1] * self.MwT:.5e}, '
                            f'{P_dpl[0, -1] * self.MwP:.5e}, {T_d[0, 0] * self.MwT:.5e}\n')

                    if self.progress_callback:
                        sample = {
                            "time": time_min,
                            "T_pl": T_pl[i] * self.MwT,
                            "PT_pl": PT_pl[i] * self.MwT,
                            "T_d": T_d[0, 0] * self.MwT,
                        }
                        self.progress_callback("sample", i, total_time_steps, sample)

                # --- Cinetica paziente ---
                alpha = self.alpha_t(Vis[i], Vpl[i])
                dVpl = -(1 - alpha) * self.Quf + self.Qr
                dVis = -alpha * self.Quf
                Vpl[i + 1] = Vpl[i] + dVpl * self.dt
                Vis[i + 1] = Vis[i] + dVis * self.dt

                dT_pl = (-self.Qp * T_pl[i] + (self.Qp - self.Quf + self.Qr) * T_tpl[0, -1] +
                         self.Kip_T * (T_is[i] - T_pl[i]) + alpha * self.Quf * T_is[i] +
                         (self.k2 * PT_pl[i] - self.k1 * P_pl[i] * T_pl[i]) * Vpl[i]) / Vpl[i] - dVpl * T_pl[i] / Vpl[i]
                T_pl[i + 1] = T_pl[i] + self.dt * dT_pl

                dPT_pl = (-self.Qp * PT_pl[i] + (self.Qp - self.Quf + self.Qr) * PT_tpl[0, -1] +
                          (self.k1 * P_pl[i] * T_pl[i] - self.k2 * PT_pl[i]) * Vpl[i]) / Vpl[i] - dVpl * PT_pl[i] / Vpl[i]
                PT_pl[i + 1] = PT_pl[i] + self.dt * dPT_pl

                dP_pl = (-self.Qp * P_pl[i] + (self.Qp - self.Quf + self.Qr) * P_tpl[0, -1] +
                         (-self.k1 * P_pl[i] * T_pl[i] + self.k2 * PT_pl[i]) * Vpl[i]) / Vpl[i] - dVpl * P_pl[i] / Vpl[i]
                P_pl[i + 1] = P_pl[i] + self.dt * dP_pl

                dT_is = (self.Kic_T * (T_ic[i] - T_is[i]) - self.Kip_T * (T_is[i] - T_pl[i]) -
                         alpha * self.Quf * T_is[i] + (self.k2 * PT_is[i] - self.k1 * P_is[i] * T_is[i]) * Vis[i]) / Vis[i] - dVis * T_is[i] / Vis[i]
                T_is[i + 1] = T_is[i] + self.dt * dT_is

                dPT_is = (self.k1 * P_is[i] * T_is[i] - self.k2 * PT_is[i]) - dVis * PT_is[i] / Vis[i]
                PT_is[i + 1] = PT_is[i] + self.dt * dPT_is

                dP_is = (-self.k1 * P_is[i] * T_is[i] + self.k2 * PT_is[i]) - dVis * P_is[i] / Vis[i]
                P_is[i + 1] = P_is[i] + self.dt * dP_is

                dT_ic = (self.G - self.Kic_T * (T_ic[i] - T_is[i])) / Vic
                T_ic[i + 1] = T_ic[i] + self.dt * dT_ic

                # --- Condizioni al contorno del dializzatore ---
                T_dpl[1, 0] = T_pl[i + 1]
                PT_dpl[1, 0] = PT_pl[i + 1]
                P_dpl[1, 0] = P_pl[i + 1]
                T_d[1, -1] = 0

                dx_Tpl = (T_dpl[0, j_dialyzer] - T_dpl[0, j_dialyzer - 1]) / self.dx
                dT_dpl = (-1 / (self.N * self.A) * (dx_Tpl * self.Qpx(j_dialyzer * self.dx) - self.Quf / self.Lfiber * T_dpl[0, j_dialyzer]) -
                          1 / (self.N * self.A * self.Lfiber) * (self.KoA * (T_dpl[0, j_dialyzer] - T_d[0, j_dialyzer])) * self.Pe / (np.exp(self.Pe) - 1) +
                          1 / (self.N * self.A * self.Lfiber) * (-(self.Quf * T_dpl[0, j_dialyzer] * self.S_T)) +
                          (-self.k1 * P_dpl[0, j_dialyzer] * T_dpl[0, j_dialyzer] + self.k2 * PT_dpl[0, j_dialyzer]))
                T_dpl[1, j_dialyzer] = T_dpl[0, j_dialyzer] + self.dt * dT_dpl

                dx_PTdpl = (PT_dpl[0, j_dialyzer] - PT_dpl[0, j_dialyzer - 1]) / self.dx
                dPT_dpl = (-1 / (self.N * self.A) * (dx_PTdpl * self.Qpx(j_dialyzer * self.dx) - self.Quf / self.Lfiber * PT_dpl[0, j_dialyzer]) +
                           (self.k1 * P_dpl[0, j_dialyzer] * T_dpl[0, j_dialyzer] - self.k2 * PT_dpl[0, j_dialyzer]))
                PT_dpl[1, j_dialyzer] = PT_dpl[0, j_dialyzer] + self.dt * dPT_dpl

                dx_Pdpl = (P_dpl[0, j_dialyzer] - P_dpl[0, j_dialyzer - 1]) / self.dx
                dP_dpl = (-1 / (self.N * self.A) * (dx_Pdpl * self.Qpx(j_dialyzer * self.dx) - self.Quf / self.Lfiber * P_dpl[0, j_dialyzer]) +
                          (-self.k1 * P_dpl[0, j_dialyzer] * T_dpl[0, j_dialyzer] + self.k2 * PT_dpl[0, j_dialyzer]))
                P_dpl[1, j_dialyzer] = P_dpl[0, j_dialyzer] + self.dt * dP_dpl

                dx_Td = (T_d[0, k_dialyzer + 1] - T_d[0, k_dialyzer - 1]) / (2 * self.dx)
                dx_Td[0] = (T_d[0, 1] - T_d[0, 0]) / self.dx
                dT_d = (1 / (self.N * self.Ad) * (dx_Td * self.Qdx(self.dx * k_dialyzer) - self.Quf / self.Lfiber * T_d[0, k_dialyzer]) +
                        1 / (self.N * self.Ad * self.Lfiber) * (self.KoA * (T_dpl[0, k_dialyzer] - T_d[0, k_dialyzer])) * self.Pe / (np.exp(self.Pe) - 1) -
                        1 / (self.N * self.Ad * self.Lfiber) * (-(self.Quf * T_dpl[0, k_dialyzer] * self.S_T)))
                T_d[1, k_dialyzer] = T_d[0, k_dialyzer] + self.dt * dT_d

                T_dpl[0, :] = T_dpl[1, :]
                PT_dpl[0, :] = PT_dpl[1, :]
                P_dpl[0, :] = P_dpl[1, :]
                T_d[0, :] = T_d[1, :]

                # --- HDF post-diluizione: tubo venoso ---
                T_tpl[1, 0] = T_dpl[1, -1] * self.Qpx(self.Lfiber) / self.Qtube
                PT_tpl[1, 0] = PT_dpl[1, -1] * self.Qpx(self.Lfiber) / self.Qtube
                P_tpl[1, 0] = P_dpl[1, -1] * self.Qpx(self.Lfiber) / self.Qtube

                dx_tTpl = (T_tpl[0, j_tube] - T_tpl[0, j_tube - 1]) / self.dx_tube
                dt_tTpl = -(self.Qtube / self.Atube) * dx_tTpl + (-self.k1 * P_tpl[0, j_tube] * T_tpl[0, j_tube] + self.k2 * PT_tpl[0, j_tube])
                T_tpl[1, j_tube] = T_tpl[0, j_tube] + self.dt * dt_tTpl

                dx_tPTpl = (PT_tpl[0, j_tube] - PT_tpl[0, j_tube - 1]) / self.dx_tube
                dt_tPTpl = -(self.Qtube / self.Atube) * dx_tPTpl + (self.k1 * P_tpl[0, j_tube] * T_tpl[0, j_tube] - self.k2 * PT_tpl[0, j_tube])
                PT_tpl[1, j_tube] = PT_tpl[0, j_tube] + self.dt * dt_tPTpl

                dx_tPpl = (P_tpl[0, j_tube] - P_tpl[0, j_tube - 1]) / self.dx_tube
                dt_tPpl = -(self.Qtube / self.Atube) * dx_tPpl + (-self.k1 * P_tpl[0, j_tube] * T_tpl[0, j_tube] + self.k2 * PT_tpl[0, j_tube])
                P_tpl[1, j_tube] = P_tpl[0, j_tube] + self.dt * dt_tPpl

                T_tpl[0, :] = T_tpl[1, :]
                PT_tpl[0, :] = PT_tpl[1, :]
                P_tpl[0, :] = P_tpl[1, :]
            else:
                if self.progress_callback:
                    self.progress_callback("done", total_time_steps, total_time_steps, None)
        finally:
            self.running = False
            if files:
                files[0].close()
                files[1].close()


# =====================================================================
# 3. INTERFACCIA GRAFICA (Tkinter)
# =====================================================================

class PbutsGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PBUTs Simulator - Cinetica tossine legate alle proteine in dialisi")
        self.geometry("1150x780")
        self.minsize(1000, 700)

        self.sim = None
        self.sim_thread = None
        self.msg_queue = queue.Queue()
        self.out_dir = os.path.abspath("./results_PBUTs_GUI")

        self._build_layout()
        self.after(150, self._poll_queue)

    # ---------------------------------------------------------
    def _build_layout(self):
        main = ttk.Frame(self, padding=10)
        main.pack(fill="both", expand=True)

        # ---- Pannello sinistro: input ----
        left = ttk.Frame(main)
        left.pack(side="left", fill="y", padx=(0, 10))

        input_frame = ttk.LabelFrame(left, text="Parametri della simulazione", padding=12)
        input_frame.pack(fill="x")

        row = 0

        ttk.Label(input_frame, text="Tempo totale di dialisi (min):").grid(row=row, column=0, sticky="w", pady=4)
        self.var_total_time = tk.StringVar(value=str(Model_Parameters_default["total_time"]))
        ttk.Entry(input_frame, textvariable=self.var_total_time, width=14).grid(row=row, column=1, pady=4)
        row += 1

        ttk.Label(input_frame, text="Tossina analizzata:").grid(row=row, column=0, sticky="w", pady=4)
        self.var_toxin = tk.StringVar(value="IS")
        toxin_combo = ttk.Combobox(input_frame, textvariable=self.var_toxin, width=12, state="readonly",
                                    values=["IS", "pCS"])
        toxin_combo.grid(row=row, column=1, pady=4)
        row += 1

        ttk.Label(input_frame, text="Filtro utilizzato:").grid(row=row, column=0, sticky="w", pady=4)
        self.var_filter = tk.StringVar(value="Fx100")
        filter_combo = ttk.Combobox(input_frame, textvariable=self.var_filter, width=12, state="readonly",
                                     values=list(Filter_models.keys()))
        filter_combo.grid(row=row, column=1, pady=4)
        row += 1

        ttk.Separator(input_frame, orient="horizontal").grid(row=row, column=0, columnspan=2, sticky="ew", pady=8)
        row += 1

        ttk.Label(input_frame, text="Portata sangue/plasmatica Qp (mL/min):").grid(row=row, column=0, sticky="w", pady=4)
        self.var_Qp = tk.StringVar(value=f"{Dialyzation_settings_default['Qp'] * 1e6:.2f}")
        ttk.Entry(input_frame, textvariable=self.var_Qp, width=14).grid(row=row, column=1, pady=4)
        row += 1

        ttk.Label(input_frame, text="Portata dialisato Qd (mL/min):").grid(row=row, column=0, sticky="w", pady=4)
        self.var_Qd = tk.StringVar(value=f"{Dialyzation_settings_default['Qd'] * 1e6:.2f}")
        ttk.Entry(input_frame, textvariable=self.var_Qd, width=14).grid(row=row, column=1, pady=4)
        row += 1

        ttk.Separator(input_frame, orient="horizontal").grid(row=row, column=0, columnspan=2, sticky="ew", pady=8)
        row += 1

        ttk.Label(input_frame, text="Concentrazione plasmatica iniziale",
                  font=("TkDefaultFont", 9, "bold")).grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1

        ttk.Label(input_frame, text="   Tossina LIBERA iniziale (mg/L):").grid(row=row, column=0, sticky="w", pady=4)
        self.var_Tfree0 = tk.StringVar(value="2.357")
        ttk.Entry(input_frame, textvariable=self.var_Tfree0, width=14).grid(row=row, column=1, pady=4)
        row += 1

        ttk.Label(input_frame, text="   Tossina TOTALE iniziale (mg/L):").grid(row=row, column=0, sticky="w", pady=4)
        self.var_Ttotal0 = tk.StringVar(value="31.02")
        ttk.Entry(input_frame, textvariable=self.var_Ttotal0, width=14).grid(row=row, column=1, pady=4)
        row += 1

        ttk.Label(input_frame, text="(totale = libera + legata alle proteine; deve essere > libera)",
                  font=("TkDefaultFont", 8), foreground="#777", wraplength=270, justify="left").grid(
            row=row, column=0, columnspan=2, sticky="w", pady=(0, 4))
        row += 1

        ttk.Separator(input_frame, orient="horizontal").grid(row=row, column=0, columnspan=2, sticky="ew", pady=8)
        row += 1

        ttk.Label(input_frame, text="Intervallo di campionamento grafico (min):").grid(row=row, column=0, sticky="w", pady=4)
        self.var_out_time = tk.StringVar(value=str(Model_Parameters_default["out_time"]))
        ttk.Entry(input_frame, textvariable=self.var_out_time, width=14).grid(row=row, column=1, pady=4)
        row += 1

        ttk.Separator(input_frame, orient="horizontal").grid(row=row, column=0, columnspan=2, sticky="ew", pady=8)
        row += 1

        adv_label = ttk.Label(input_frame, text="Parametri avanzati (opzionali)", font=("TkDefaultFont", 9, "italic"))
        adv_label.grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1

        ttk.Label(input_frame, text="Ultrafiltrazione Quf (mL/min):").grid(row=row, column=0, sticky="w", pady=4)
        self.var_Quf = tk.StringVar(value=f"{Dialyzation_settings_default['Quf'] * 1e6:.2f}")
        ttk.Entry(input_frame, textvariable=self.var_Quf, width=14).grid(row=row, column=1, pady=4)
        row += 1

        ttk.Label(input_frame, text="Reinfusione Qr (mL/min):").grid(row=row, column=0, sticky="w", pady=4)
        self.var_Qr = tk.StringVar(value=f"{Dialyzation_settings_default['Qr'] * 1e6:.2f}")
        ttk.Entry(input_frame, textvariable=self.var_Qr, width=14).grid(row=row, column=1, pady=4)
        row += 1

        ttk.Label(input_frame, text="Passo temporale dt (min):").grid(row=row, column=0, sticky="w", pady=4)
        self.var_dt = tk.StringVar(value=str(Model_Parameters_default["dt"]))
        ttk.Entry(input_frame, textvariable=self.var_dt, width=14).grid(row=row, column=1, pady=4)
        row += 1

        self.var_save_file = tk.BooleanVar(value=True)
        ttk.Checkbutton(input_frame, text="Salva risultati su file (.txt)", variable=self.var_save_file).grid(
            row=row, column=0, columnspan=2, sticky="w", pady=(6, 0))
        row += 1

        # ---- Bottoni ----
        btn_frame = ttk.Frame(left)
        btn_frame.pack(fill="x", pady=12)
        self.btn_start = ttk.Button(btn_frame, text="Avvia simulazione", command=self._on_start)
        self.btn_start.pack(side="left", fill="x", expand=True, padx=(0, 4))
        self.btn_stop = ttk.Button(btn_frame, text="Interrompi", command=self._on_stop, state="disabled")
        self.btn_stop.pack(side="left", fill="x", expand=True, padx=(4, 0))

        self.progress = ttk.Progressbar(left, orient="horizontal", mode="determinate", length=280)
        self.progress.pack(fill="x", pady=(4, 4))

        self.status_var = tk.StringVar(value="Pronto.")
        ttk.Label(left, textvariable=self.status_var, wraplength=280, foreground="#333").pack(fill="x", pady=(2, 8))

        info = ("Nota: con dt molto piccolo (es. 1e-5 min) e tempi lunghi il numero di passi "
                "puo' essere elevato e la simulazione richiedere diversi minuti.")
        ttk.Label(left, text=info, wraplength=280, foreground="#777",
                  font=("TkDefaultFont", 8)).pack(fill="x")

        # ---- Pannello destro: grafici ----
        right = ttk.Frame(main)
        right.pack(side="left", fill="both", expand=True)

        self.fig = Figure(figsize=(7.5, 7), dpi=100)
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        self._reset_axes()

        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def _reset_axes(self):
        self.ax1.clear()
        self.ax1.set_title("Concentrazioni plasmatiche (aggiornamento in tempo reale)")
        self.ax1.set_xlabel("Tempo (min)")
        self.ax1.set_ylabel("Concentrazione (mg/L)")
        self.ax1.grid(True, alpha=0.3)
        (self.line_free,) = self.ax1.plot([], [], marker="o", markersize=4, linewidth=2,
                                           color="#1f77b4", label="Tossina LIBERA")
        (self.line_total,) = self.ax1.plot([], [], marker="s", markersize=4, linewidth=2,
                                            color="#2ca02c", label="Tossina TOTALE (libera+legata)")
        (self.line_bound,) = self.ax1.plot([], [], marker="", linewidth=1, linestyle=":",
                                            color="#ff7f0e", alpha=0.6, label="Legata alle proteine")
        self.ax1.legend(loc="best", fontsize=8)

        self.ax2.clear()
        self.ax2.set_title("Tossina libera nel dialisato in uscita")
        self.ax2.set_xlabel("Tempo (min)")
        self.ax2.set_ylabel("Concentrazione (mg/L)")
        self.ax2.grid(True, alpha=0.3)
        (self.line_dial,) = self.ax2.plot([], [], marker="o", markersize=4, linewidth=2, color="#9467bd")

        self.fig.tight_layout()
        if hasattr(self, "canvas"):
            self.canvas.draw()

        # buffer dati per l'aggiornamento incrementale
        self._plot_t = []
        self._plot_free = []
        self._plot_total = []
        self._plot_bound = []
        self._plot_dial = []

    # ---------------------------------------------------------
    def _read_inputs(self):
        try:
            total_time = float(self.var_total_time.get())
            dt = float(self.var_dt.get())
            out_time = float(self.var_out_time.get())
            Qp = float(self.var_Qp.get()) * 1e-6     # mL/min -> m^3/min
            Qd = float(self.var_Qd.get()) * 1e-6
            Quf = float(self.var_Quf.get()) * 1e-6
            Qr = float(self.var_Qr.get()) * 1e-6
            T_free0 = float(self.var_Tfree0.get())    # mg/L
            T_total0 = float(self.var_Ttotal0.get())  # mg/L
            toxin = self.var_toxin.get()
            filt = self.var_filter.get()
        except ValueError:
            raise ValueError("Controlla i valori numerici inseriti (usa il punto come separatore decimale).")

        if total_time <= 0 or dt <= 0:
            raise ValueError("Tempo totale e passo temporale devono essere numeri positivi.")
        if out_time <= 0:
            raise ValueError("L'intervallo di campionamento del grafico deve essere positivo.")
        if out_time > total_time:
            raise ValueError("L'intervallo di campionamento non puo' essere maggiore del tempo totale di dialisi.")
        if Qp <= 0 or Qd <= 0:
            raise ValueError("Le portate di sangue e dialisato devono essere positive.")
        if toxin not in Toxin_parameters:
            raise ValueError("Tossina non valida.")
        if filt not in Filter_models:
            raise ValueError("Filtro non valido.")
        if T_free0 < 0 or T_total0 < 0:
            raise ValueError("Le concentrazioni iniziali di tossina non possono essere negative.")
        if T_total0 <= T_free0:
            raise ValueError("La concentrazione TOTALE iniziale deve essere maggiore di quella LIBERA "
                              "(la differenza rappresenta la quota legata alle proteine).")

        n_steps = int(total_time / dt)
        if n_steps > 30_000_000:
            raise ValueError(
                f"La combinazione tempo/dt richiederebbe {n_steps:,} passi: troppi per una simulazione "
                f"interattiva. Aumenta dt o riduci il tempo totale.")

        model_params = dict(Model_Parameters_default)
        model_params["total_time"] = total_time
        model_params["dt"] = dt
        model_params["out_time"] = out_time

        dialysis_settings = {"Qp": Qp, "Quf": Quf, "Qd": Qd, "Qr": Qr}

        # Costruzione parametri paziente a partire dalle concentrazioni iniziali
        # scelte dall'utente (libera e totale), convertite in mol/m^3.
        MwT = Toxin_parameters[toxin]["MwT"]
        Tpl0 = T_free0 / MwT
        PTpl0 = (T_total0 - T_free0) / MwT
        patient = dict(Patient_default)
        patient["Tpl0"] = Tpl0
        patient["PTpl0"] = PTpl0
        patient["Tis0"] = Tpl0        # tossina libera: stessa concentrazione iniziale di plasma e interstizio
        patient["PTis0"] = PTpl0 / 2  # legata: meta' della concentrazione plasmatica (come nel modello originale)
        patient["Tic0"] = Tpl0

        return toxin, filt, model_params, dialysis_settings, patient, n_steps

    # ---------------------------------------------------------
    def _on_start(self):
        if self.sim_thread and self.sim_thread.is_alive():
            return
        try:
            toxin, filt, model_params, dialysis_settings, patient, n_steps = self._read_inputs()
        except ValueError as e:
            messagebox.showerror("Parametri non validi", str(e))
            return

        self._reset_axes()
        self.progress["value"] = 0
        self.progress["maximum"] = n_steps
        self.status_var.set(f"Avvio simulazione ({n_steps:,} passi)...")
        self.btn_start.config(state="disabled")
        self.btn_stop.config(state="normal")

        self.sim = PbutsPredictor(
            filter_model=filt,
            Toxin=toxin,
            Patient=patient,
            Model_Parameters=model_params,
            Dialyzation_settings=dialysis_settings,
            out_name="GUI",
            out_path=self.out_dir,
            progress_callback=self._progress_callback,
            save_to_file=self.var_save_file.get(),
        )

        for msg in self.sim.stability_msgs:
            self.msg_queue.put(("warning", msg))

        self.sim_thread = threading.Thread(target=self.sim.run, daemon=True)
        self.sim_thread.start()

    def _on_stop(self):
        if self.sim:
            self.sim.stop()
        self.status_var.set("Interruzione in corso...")
        self.btn_stop.config(state="disabled")

    def _progress_callback(self, kind, i, total, message):
        # Chiamato dal thread di simulazione: passiamo tutto alla coda
        # thread-safe, la GUI la legge nel main thread con _poll_queue.
        self.msg_queue.put((kind, (i, total, message)))

    # ---------------------------------------------------------
    def _poll_queue(self):
        try:
            while True:
                kind, payload = self.msg_queue.get_nowait()
                if kind == "warning":
                    self.status_var.set(payload)
                elif kind == "sample":
                    i, total, data = payload
                    self.progress["value"] = i
                    self.status_var.set(
                        f"Simulazione in corso: t={data['time']:.1f} min "
                        f"({100 * i / total:.1f}%) - libera={data['T_pl']:.3f} mg/L, "
                        f"totale={data['T_pl'] + data['PT_pl']:.3f} mg/L")
                    self._append_point(data)
                elif kind in ("done", "stopped", "error"):
                    i, total, message = payload
                    self.progress["value"] = total if kind == "done" else i
                    self.btn_start.config(state="normal")
                    self.btn_stop.config(state="disabled")
                    if kind == "done":
                        self.status_var.set("Simulazione completata.")
                    elif kind == "stopped":
                        self.status_var.set(message or "Simulazione interrotta.")
                    else:
                        self.status_var.set(message or "Errore durante la simulazione.")
                        messagebox.showerror("Errore simulazione", message or "Errore sconosciuto.")
        except queue.Empty:
            pass
        self.after(150, self._poll_queue)

    # ---------------------------------------------------------
    def _append_point(self, data):
        """Aggiunge un nuovo punto campionato e ridisegna le curve in tempo reale."""
        t = data["time"]
        free = data["T_pl"]
        bound = data["PT_pl"]
        total = free + bound
        dial = data["T_d"]

        self._plot_t.append(t)
        self._plot_free.append(free)
        self._plot_bound.append(bound)
        self._plot_total.append(total)
        self._plot_dial.append(dial)

        title = "Concentrazioni plasmatiche"
        if self.sim is not None:
            title += f" - {self.sim.toxin_name} - filtro {self.sim.filter_model}"
        self.ax1.set_title(title)

        self.line_free.set_data(self._plot_t, self._plot_free)
        self.line_total.set_data(self._plot_t, self._plot_total)
        self.line_bound.set_data(self._plot_t, self._plot_bound)
        self.ax1.relim()
        self.ax1.autoscale_view()

        self.line_dial.set_data(self._plot_t, self._plot_dial)
        self.ax2.relim()
        self.ax2.autoscale_view()

        self.fig.tight_layout()
        self.canvas.draw_idle()


# =====================================================================
# 4. AVVIO APPLICAZIONE
# =====================================================================

if __name__ == "__main__":
    app = PbutsGUI()
    app.mainloop()
