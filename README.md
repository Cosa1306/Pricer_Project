Option Pricer with Streamlit
This project is an Option Pricer that allows you to calculate the price of an option based on different financial models. Originally developed with Tkinter, the pricer has now been migrated to Streamlit for a more modern and interactive web interface. The pricer now supports several option pricing models:

Black-Scholes: The classic model for pricing European options.

Binomial Tree: A discrete approach for pricing American options.

Monte Carlo Simulation: A simulation-based method used for pricing Asian options.
The tool allows you to choose between these different methods, specify option parameters (such as stock price, volatility, risk-free rate, etc.), and display the results in real-time via an interactive graphical interface on Streamlit. Additionally, it offers the calculation of Greeks (Delta, Gamma, Vega, Theta, Rho) to analyze option price sensitivities.

In addition to the option pricer, we have also designed a bond pricer that allows users to calculate essential bond metrics such as the clean price, dirty price, and accrued interest, as well as advanced analytics like the Macauley duration, modified duration, and convexity for assessing interest rate sensitivity and price risk.
