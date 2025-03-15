```python
"""
@author: Sonny Genovese & Thomas-Victor Coll
"""
import streamlit as st
from datetime import datetime, timedelta

def calculate_dirty_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, yield_to_maturity):
    coupon_dates = []
    next_coupon_date = last_coupon_date + timedelta(days=365 / frequency)
    while next_coupon_date <= maturity_date:
        if next_coupon_date > today_date:
            coupon_dates.append(next_coupon_date)
        next_coupon_date += timedelta(days=365 / frequency)

    present_value = 0
    for date in coupon_dates:
        time_to_payment = (date - today_date).days / 365.0
        discount_factor = (1 + yield_to_maturity / frequency) ** (frequency * time_to_payment)
        present_value += (coupon_rate * face_value / frequency) / discount_factor

    time_to_maturity = (maturity_date - today_date).days / 365.0
    discount_factor = (1 + yield_to_maturity / frequency) ** (frequency * time_to_maturity)
    present_value += face_value / discount_factor

    return present_value

def calculate_accrued_interest(coupon_rate, frequency, last_coupon_date, today_date, face_value):
    days_since_last_coupon = (today_date - last_coupon_date).days
    next_coupon_date = last_coupon_date + timedelta(days=365 / frequency)
    days_in_period = (next_coupon_date - last_coupon_date).days
    accrued_interest = (coupon_rate * face_value / frequency) * (days_since_last_coupon / days_in_period)
    return accrued_interest

def calculate_clean_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, yield_to_maturity):
    dirty_price = calculate_dirty_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, yield_to_maturity)
    accrued_interest = calculate_accrued_interest(coupon_rate, frequency, last_coupon_date, today_date, face_value)
    clean_price = dirty_price - accrued_interest
    return clean_price

# Streamlit UI
st.title("Bond Pricer")

coupon_rate = st.number_input("Coupon Rate (%)", min_value=0.0, max_value=100.0, value=5.0, step=0.1) / 100
ytm = st.number_input("Yield to Maturity (YTM) (%)", min_value=0.0, max_value=100.0, value=3.0, step=0.1) / 100
frequency = st.number_input("Frequency (Payments per Year)", min_value=1, max_value=12, value=2, step=1)
last_coupon_date = st.date_input("Last Coupon Date", datetime.today() - timedelta(days=180))
today_date = st.date_input("Today's Date", datetime.today())
maturity_date = st.date_input("Maturity Date", datetime.today() + timedelta(days=365 * 10))
face_value = 1000  # Fixed face value

if st.button("Calculate Bond Prices"):
    clean_price = calculate_clean_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, ytm)
    dirty_price = calculate_dirty_price(coupon_rate, frequency, last_coupon_date, today_date, maturity_date, face_value, ytm)
    accrued_interest = calculate_accrued_interest(coupon_rate, frequency, last_coupon_date, today_date, face_value)
    
    st.write(f"**Clean Price:** {clean_price:.4f}")
    st.write(f"**Dirty Price:** {dirty_price:.4f}")
    st.write(f"**Accrued Interest:** {accrued_interest:.4f}")

tk.Label(root, textvariable=results, justify="left").grid(row=7, column=0, columnspan=2)

# Start application
root.mainloop()
