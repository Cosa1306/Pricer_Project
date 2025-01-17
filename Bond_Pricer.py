```python
"""
@author: Sonny Genovese & Thomas-Victor Coll
"""
import tkinter as tk
from tkinter.messagebox import showinfo
from datetime import datetime, timedelta

def calculate_dirty_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, yield_to_maturity):
    # Generate future coupon payment dates (only future)
    coupon_dates = []
    next_coupon_date = last_coupon_date + timedelta(days=365 / frequency)
    while next_coupon_date <= maturity_date:
        if next_coupon_date > today_date:
            coupon_dates.append(next_coupon_date)
        next_coupon_date += timedelta(days=365 / frequency)

    # Calculate present value of all future cash flows
    present_value = 0
    for date in coupon_dates:
        time_to_payment = (date - today_date).days / 365.0
        discount_factor = (1 + yield_to_maturity / frequency) ** (frequency * time_to_payment)
        present_value += (coupon_rate * face_value / frequency) / discount_factor

    # Add discounted face value at maturity
    time_to_maturity = (maturity_date - today_date).days / 365.0
    discount_factor = (1 + yield_to_maturity / frequency) ** (frequency * time_to_maturity)
    present_value += face_value / discount_factor

    return present_value

def calculate_accrued_interest(coupon_rate, frequency, last_coupon_date, today_date, face_value):
    # Calculate days since last coupon and days in period
    days_since_last_coupon = (today_date - last_coupon_date).days
    next_coupon_date = last_coupon_date + timedelta(days=365 / frequency)
    days_in_period = (next_coupon_date - last_coupon_date).days
    accrued_interest = (coupon_rate * face_value / frequency) * (days_since_last_coupon / days_in_period)
    return accrued_interest

def calculate_clean_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, yield_to_maturity):
    # Calculate dirty price and subtract accrued interest to get clean price
    dirty_price = calculate_dirty_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, yield_to_maturity)
    accrued_interest = calculate_accrued_interest(coupon_rate, frequency, last_coupon_date, today_date, face_value)
    clean_price = dirty_price - accrued_interest
    return clean_price

def calculate_and_display():
    try:
        # Fetch inputs from GUI
        coupon_rate = float(coupon_rate_var.get()) / 100  # Convert % to decimal
        yield_to_maturity = float(ytm_var.get()) / 100  # Convert % to decimal
        frequency = int(frequency_var.get())
        last_coupon_date = datetime.strptime(last_coupon_var.get(), "%Y-%m-%d")
        today_date = datetime.strptime(today_date_var.get(), "%Y-%m-%d")
        maturity_date = datetime.strptime(maturity_date_var.get(), "%Y-%m-%d")
        face_value = 1000  # Assuming a standard face value of 1000

        # Calculate prices
        clean_price = calculate_clean_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, yield_to_maturity)
        dirty_price = calculate_dirty_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, yield_to_maturity)
        accrued_interest = calculate_accrued_interest(coupon_rate, frequency, last_coupon_date, today_date, face_value)

        # Update results in GUI
        results.set(
            f"Clean Price: {clean_price:.4f}\n"
            f"Dirty Price: {dirty_price:.4f}\n"
            f"Accrued Interest: {accrued_interest:.4f}"
        )
    except Exception as e:
        showinfo("Error", f"An error occurred: {e}")

# GUI Setup
root = tk.Tk()
root.title("Bond Pricer")

coupon_rate_var = tk.StringVar()
ytm_var = tk.StringVar()
frequency_var = tk.StringVar()
last_coupon_var = tk.StringVar()
today_date_var = tk.StringVar()
maturity_date_var = tk.StringVar()
results = tk.StringVar(value="")

tk.Label(root, text="Coupon Rate (%):").grid(row=0, column=0, sticky="e")
tk.Entry(root, textvariable=coupon_rate_var).grid(row=0, column=1)

tk.Label(root, text="YTM (%):").grid(row=1, column=0, sticky="e")
tk.Entry(root, textvariable=ytm_var).grid(row=1, column=1)

tk.Label(root, text="Frequency :").grid(row=2, column=0, sticky="e")
tk.Entry(root, textvariable=frequency_var).grid(row=2, column=1)

tk.Label(root, text="Last Coupon Date (YYYY-MM-DD):").grid(row=3, column=0, sticky="e")
tk.Entry(root, textvariable=last_coupon_var).grid(row=3, column=1)

tk.Label(root, text="Today Date (YYYY-MM-DD):").grid(row=4, column=0, sticky="e")
tk.Entry(root, textvariable=today_date_var).grid(row=4, column=1)

tk.Label(root, text="Maturity Date (YYYY-MM-DD):").grid(row=5, column=0, sticky="e")
tk.Entry(root, textvariable=maturity_date_var).grid(row=5, column=1)

tk.Button(root, text="Calculate", command=calculate_and_display).grid(row=6, column=0, columnspan=2)

tk.Label(root, textvariable=results, justify="left").grid(row=7, column=0, columnspan=2)

# Start application
root.mainloop()
