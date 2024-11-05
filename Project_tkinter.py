#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 09:35:32 2024

@author: Sonny Genovese et Thomas-Victor Coll
"""

from tkinter import *
import numpy as np
from scipy.stats import norm


# Definition formule Black & Scholes
def call_pricer(stock_price, strike_price, r_free, maturity, sigma, option_name):
    d1 = (np.log(stock_price / strike_price) + (r_free + (sigma**2) / 2) * maturity) / (sigma * np.sqrt(maturity))
    d2 = d1 - sigma * np.sqrt(maturity)

    N1 = norm.cdf(d1)
    N2 = norm.cdf(d2)
    price = stock_price * N1 - strike_price * np.exp(-r_free * maturity) * N2
    return price


def put_pricer(stock_price, strike_price, r_free, maturity, sigma, option_name):
    d1 = (np.log(stock_price / strike_price) + (r_free + (sigma**2) / 2) * maturity) / (sigma * np.sqrt(maturity))
    d2 = d1 - sigma * np.sqrt(maturity)

    N1 = norm.cdf(d1)
    N2 = norm.cdf(d2)
    price = price = strike_price * np.exp(-r_free * maturity) * (1 - N2) - stock_price * (1 - N1)
    return price

# Definition bond pricer
def bond_pricer(m, t, ytm, bond_fv, coupon_r):
  bond_price=((bond_fv*coupon_r/m*(1-(1+ytm/m)**(-m*t)))/(ytm/m))+bond_fv*(1+(ytm/m))**(-m*t)
  return bond_price

# ratios
def sharp_ratio(r,r_free,sigma):
  sr=(r-r_free)/sigma
  return sr

# duration
def duration(d,c,m,y,n,t):
  numerator = sum([(t*c)/((1+y)**t)+((n*m)/(1+y)**n)])
  denominator = sum([(c/(1+y)**t)+m/(1+y)**n])
  return numerator/denominator

#Main window
fenetre = Tk()
fenetre.title("Pricer")

label = Label(fenetre, text="Bonjour! Que souhaitez-vous calculer?")
label.pack()

# List of tools
liste = Listbox(fenetre)
liste.insert(1, "option call")
liste.insert(2, "option put")
liste.insert(3, "obligation")
liste.insert(4, "sharp_ratio")
liste.insert(5, "duration")
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
    selected_index = liste.curselection()
    if selected_index:
        selected_item = liste.get(selected_index)

        if selected_item == "option call":
            call_window = Toplevel(fenetre)
            call_window.title("Option Call Pricing Parameters")
            entry_S = create_labeled_entry(call_window, "Stock price (S):")
            entry_K = create_labeled_entry(call_window, "Strike price (K):")
            entry_r = create_labeled_entry(call_window, "Risk-free rate (r):")
            entry_T = create_labeled_entry(call_window, "Time to maturity (T):")
            entry_sigma = create_labeled_entry(call_window, "Volatility (sigma):")

            def calculate_call_price():
                S = float(entry_S.get())
                K = float(entry_K.get())
                r = float(entry_r.get())
                T = float(entry_T.get())
                sigma = float(entry_sigma.get())
                price = call_pricer(S, K, r, T, sigma)
                result_label.config(text=f"Option Call Price: {price:.2f}")

            calculate_button = Button(call_window, text="Calculate", command=calculate_call_price)
            calculate_button.pack()
            result_label = Label(call_window, text="")
            result_label.pack()

        elif selected_item == "option put":
            put_window = Toplevel(fenetre)
            put_window.title("Option Put Pricing Parameters")
            entry_S = create_labeled_entry(put_window, "Stock price (S):")
            entry_K = create_labeled_entry(put_window, "Strike price (K):")
            entry_r = create_labeled_entry(put_window, "Risk-free rate (r):")
            entry_T = create_labeled_entry(put_window, "Time to maturity (T):")
            entry_sigma = create_labeled_entry(put_window, "Volatility (sigma):")

            def calculate_put_price():
                S = float(entry_S.get())
                K = float(entry_K.get())
                r = float(entry_r.get())
                T = float(entry_T.get())
                sigma = float(entry_sigma.get())
                price = put_pricer(S, K, r, T, sigma)
                result_label.config(text=f"Option Put Price: {price:.2f}")

            calculate_button = Button(put_window, text="Calculate", command=calculate_put_price)
            calculate_button.pack()
            result_label = Label(put_window, text="")
            result_label.pack()

        elif selected_item == "obligation":
            bond_window = Toplevel(fenetre)
            bond_window.title("Bond Pricing Parameters")
            entry_m = create_labeled_entry(bond_window, "Number of payments per year (m):")
            entry_t = create_labeled_entry(bond_window, "Total years to maturity (t):")
            entry_ytm = create_labeled_entry(bond_window, "Yield to maturity (ytm):")
            entry_bond_fv = create_labeled_entry(bond_window, "Bond face value (FV):")
            entry_coupon_r = create_labeled_entry(bond_window, "Coupon rate (coupon_r):")

            def calculate_bond_price():
                m = float(entry_m.get())
                t = float(entry_t.get())
                ytm = float(entry_ytm.get())
                bond_fv = float(entry_bond_fv.get())
                coupon_r = float(entry_coupon_r.get())
                price = bond_pricer(m, t, ytm, bond_fv, coupon_r)
                result_label.config(text=f"Bond Price: {price:.2f}")

            calculate_button = Button(bond_window, text="Calculate", command=calculate_bond_price)
            calculate_button.pack()
            result_label = Label(bond_window, text="")
            result_label.pack()

        elif selected_item == "sharp_ratio":
            bond_window = Toplevel(fenetre)
            bond_window.title("Sharp ratio calculator")
            entry_m = create_labeled_entry(bond_window, "Number of payments per year (m):")
            entry_t = create_labeled_entry(bond_window, "Total years to maturity (t):")
            entry_ytm = create_labeled_entry(bond_window, "Yield to maturity (ytm):")
            entry_bond_fv = create_labeled_entry(bond_window, "Bond face value (FV):")
            entry_coupon_r = create_labeled_entry(bond_window, "Coupon rate (coupon_r):")


# Bind the selection event
liste.bind("<<ListboxSelect>>", on_select)

fenetre.mainloop()
