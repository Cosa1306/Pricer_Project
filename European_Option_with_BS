#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 09:35:32 2024

@author: Sonny Genovese & Thomas-Victor Coll
"""

from tkinter import *
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Black & Scholes formula definition
def call_pricer(stock_price, strike_price, r_free, maturity, sigma):
    d1 = (np.log(stock_price / strike_price) + (r_free + (sigma**2) / 2) * maturity) / (sigma * np.sqrt(maturity))
    d2 = d1 - sigma * np.sqrt(maturity)
    N1 = norm.cdf(d1)
    N2 = norm.cdf(d2)
    price = stock_price * N1 - strike_price * np.exp(-r_free * maturity) * N2
    return price

def put_pricer(stock_price, strike_price, r_free, maturity, sigma):
    d1 = (np.log(stock_price / strike_price) + (r_free + (sigma**2) / 2) * maturity) / (sigma * np.sqrt(maturity))
    d2 = d1 - sigma * np.sqrt(maturity)
    N1 = norm.cdf(d1)
    N2 = norm.cdf(d2)
    price = strike_price * np.exp(-r_free * maturity) * (1 - N2) - stock_price * (1 - N1)
    return price

# Main window
fenetre = Tk()
fenetre.title("Pricer")
fenetre.geometry("400x300")

label = Label(fenetre, text="Hello! What would you like to calculate?")
label.pack()

# List of tools
liste = Listbox(fenetre)
liste.insert(1, "option call")
liste.insert(2, "option put")
liste.pack()

# Helper function to create labeled entry
def create_labeled_entry(window, text):
    label = Label(window, text=text)
    label.pack()
    entry = Entry(window)
    entry.pack()
    return entry

# Function to handle selection and open the respective input window
def on_select(event):
    selection = event.widget.get(event.widget.curselection())
    if selection == "option call":
        create_option_input_window("call")
    elif selection == "option put":
        create_option_input_window("put")

def create_option_input_window(option_type):
    window = Toplevel(fenetre)
    window.geometry("400x300")
    window.title(f"{option_type.capitalize()} Option Pricing Parameters")
    entry_S = create_labeled_entry(window, "Stock price (S):")
    entry_K = create_labeled_entry(window, "Strike price (K):")
    entry_r = create_labeled_entry(window, "Risk-free rate (r):")
    entry_T = create_labeled_entry(window, "Time to maturity (T):")
    entry_sigma = create_labeled_entry(window, "Volatility (sigma):")
    result_label = Label(window, text="")
    result_label.pack()

    # Checkboxes for Greeks
    delta_var = BooleanVar()
    gamma_var = BooleanVar()
    vega_var = BooleanVar()
    theta_var = BooleanVar()
    rho_var = BooleanVar()

    create_checkbox(window, "Delta", delta_var)
    create_checkbox(window, "Gamma", gamma_var)
    create_checkbox(window, "Vega", vega_var)
    create_checkbox(window, "Theta", theta_var)
    create_checkbox(window, "Rho", rho_var)

    if option_type == "call":
        calculate_button = Button(window, text="Calculate", command=lambda: calculate_price(entry_S, entry_K, entry_r, entry_T, entry_sigma, result_label, call_pricer, delta_var, gamma_var, vega_var, theta_var, rho_var, window))
    else:
        calculate_button = Button(window, text="Calculate", command=lambda: calculate_price(entry_S, entry_K, entry_r, entry_T, entry_sigma, result_label, put_pricer, delta_var, gamma_var, vega_var, theta_var, rho_var, window))
    calculate_button.pack()

def calculate_price(entry_S, entry_K, entry_r, entry_T, entry_sigma, result_label, pricer_func, delta_var, gamma_var, vega_var, theta_var, rho_var, window):
    try:
        S = float(entry_S.get())
        K = float(entry_K.get())
        r = float(entry_r.get())
        T = float(entry_T.get())
        sigma = float(entry_sigma.get())
        if S <= 0 or K <= 0 or r < 0 or T <= 0 or sigma <= 0:
            raise ValueError("Input values must be positive numbers")
        price = pricer_func(S, K, r, T, sigma)
        result_label.config(text=f"Option Price: {price:.2f}")

        # Calculate Greeks if checkboxes are selected
        if delta_var.get():
            delta = calculate_delta(S, K, r, T, sigma, pricer_func)
            result_label.config(text=result_label.cget("text") + f"\nDelta: {delta:.4f}")
        if gamma_var.get():
            gamma = calculate_gamma(S, K, r, T, sigma)
            result_label.config(text=result_label.cget("text") + f"\nGamma: {gamma:.4f}")
        if vega_var.get():
            vega = calculate_vega(S, K, r, T, sigma)
            result_label.config(text=result_label.cget("text") + f"\nVega: {vega:.4f}")
        if theta_var.get():
            theta = calculate_theta(S, K, r, T, sigma, pricer_func)
            result_label.config(text=result_label.cget("text") + f"\nTheta: {theta:.4f}")
        if rho_var.get():
            rho = calculate_rho(S, K, r, T, sigma, pricer_func)
            result_label.config(text=result_label.cget("text") + f"\nRho: {rho:.4f}")

        # Add a button to show the plot if it doesn't already exist
        if not hasattr(window, 'plot_button'):
            plot_button = Button(window, text="Show Greeks Plot", command=lambda: show_greeks_plot(S, K, r, T, sigma, pricer_func, delta_var, gamma_var, vega_var, theta_var, rho_var))
            plot_button.pack()
            window.plot_button = plot_button  # Store the button in the window object

    except ValueError as e:
        result_label.config(text=str(e))

def create_checkbox(window, text, var):
    checkbox = Checkbutton(window, text=text, variable=var)
    checkbox.pack()

def calculate_delta(S, K, r, T, sigma, pricer_func):
    d1 = (np.log(S / K) + (r + (sigma**2) / 2) * T) / (sigma * np.sqrt(T))
    if pricer_func == call_pricer:
        return norm.cdf(d1)
    else:
        return norm.cdf(d1) - 1

def calculate_gamma(S, K, r, T, sigma):
    d1 = (np.log(S / K) + (r + (sigma**2) / 2) * T) / (sigma * np.sqrt(T))
    return norm.pdf(d1) / (S * sigma * np.sqrt(T))

def calculate_vega(S, K, r, T, sigma):
    d1 = (np.log(S / K) + (r + (sigma**2) / 2) * T) / (sigma * np.sqrt(T))
    return S * norm.pdf(d1) * np.sqrt(T)

def calculate_theta(S, K, r, T, sigma, pricer_func):
    d1 = (np.log(S / K) + (r + (sigma**2) / 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    if pricer_func == call_pricer:
        return -(S * norm.pdf(d1) * sigma) / (2 * np.sqrt(T)) - r * K * np.exp(-r * T) * norm.cdf(d2)
    else:
        return -(S * norm.pdf(d1) * sigma) / (2 * np.sqrt(T)) + r * K * np.exp(-r * T) * norm.cdf(-d2)

def calculate_rho(S, K, r, T, sigma, pricer_func):
    d1 = (np.log(S / K) + (r + (sigma**2) / 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    if pricer_func == call_pricer:
        return K * T * np.exp(-r * T) * norm.cdf(d2)
    else:
        return -K * T * np.exp(-r * T) * norm.cdf(-d2)

def show_greeks_plot(S, K, r, T, sigma, pricer_func, delta_var, gamma_var, vega_var, theta_var, rho_var):
    plot_window = Toplevel(fenetre)
    plot_window.title("Greeks Plot")

    canvas = Canvas(plot_window)
    scrollbar = Scrollbar(plot_window, orient="vertical", command=canvas.yview)
    scrollable_frame = Frame(canvas)

    scrollable_frame.bind(
        "<Configure>",
        lambda e: canvas.configure(
            scrollregion=canvas.bbox("all")
        )
    )

    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)

    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

    row, col = 0, 0

    if delta_var.get():
        stock_prices = np.linspace(S * 0.5, S * 1.5, 100)
        delta_values = [calculate_delta(S, K, r, T, sigma, pricer_func) for S in stock_prices]
        fig, ax = plt.subplots()
        ax.plot(stock_prices, delta_values, label='Delta', color='blue')
        ax.set_xlabel('Price of the underlying asset (S)', fontsize=14)
        ax.set_ylabel('Delta (∂V/∂S)', fontsize=14)
        ax.legend(fontsize=14)
        ax.grid(True)
        plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.2)
        canvas = FigureCanvasTkAgg(fig, master=scrollable_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=row, column=col, padx=5, pady=5)
        col += 1
        if col > 2:
            col = 0
            row += 1

    if gamma_var.get():
        stock_prices = np.linspace(S * 0.5, S * 1.5, 100)
        gamma_values = [calculate_gamma(S, K, r, T, sigma) for S in stock_prices]
        fig, ax = plt.subplots()
        ax.plot(stock_prices, gamma_values, label='Gamma', color='green')
        ax.set_xlabel('Price of the underlying asset (S)', fontsize=14)
        ax.set_ylabel('Gamma (∂²V/∂S²)', fontsize=14)
        ax.legend(fontsize=14)
        ax.grid(True)
        plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.2)
        canvas = FigureCanvasTkAgg(fig, master=scrollable_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=row, column=col, padx=5, pady=5)
        col += 1
        if col > 2:
            col = 0
            row += 1

    if vega_var.get():
        volatilities = np.linspace(sigma * 0.5, sigma * 1.5, 100)
        vega_values = [calculate_vega(S, K, r, T, sigma) for sigma in volatilities]
        fig, ax = plt.subplots()
        ax.plot(volatilities, vega_values, label='Vega', color='red')
        ax.set_xlabel('Implied volatility (σ)', fontsize=14)
        ax.set_ylabel('Vega (∂V/∂σ)', fontsize=14)
        ax.legend(fontsize=14)
        ax.grid(True)
        plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.2)
        canvas = FigureCanvasTkAgg(fig, master=scrollable_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=row, column=col, padx=5, pady=5)
        col += 1
        if col > 2:
            col = 0
            row += 1

    if theta_var.get():
        times = np.linspace(T * 0.5, T * 1.5, 100)
        theta_values = [calculate_theta(S, K, r, T, sigma, pricer_func) for T in times]
        fig, ax = plt.subplots()
        ax.plot(times, theta_values, label='Theta', color='purple')
        ax.set_xlabel('Time to expiration (T)', fontsize=14)
        ax.set_ylabel('Theta (∂V/∂T)', fontsize=14)
        ax.legend(fontsize=14)
        ax.grid(True)
        plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.2)
        canvas = FigureCanvasTkAgg(fig, master=scrollable_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=row, column=col, padx=5, pady=5)
        col += 1
        if col > 2:
            col = 0
            row += 1

    if rho_var.get():
        interest_rates = np.linspace(r * 0.5, r * 1.5, 100)
        rho_values = [calculate_rho(S, K, r, T, sigma, pricer_func) for r in interest_rates]
        fig, ax = plt.subplots()
        ax.plot(interest_rates, rho_values, label='Rho', color='orange')
        ax.set_xlabel('Interest rate (r)', fontsize=14)
        ax.set_ylabel('Rho (∂V/∂r)', fontsize=14)
        ax.legend(fontsize=14)
        ax.grid(True)
        plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.2)
        canvas = FigureCanvasTkAgg(fig, master=scrollable_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=row, column=col, padx=5, pady=5)
        col += 1
        if col > 2:
            col = 0
            row += 1

# Bind the selection event
liste.bind("<<ListboxSelect>>", on_select)

fenetre.mainloop()
