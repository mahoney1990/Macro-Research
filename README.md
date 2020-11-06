# Macro-Research
The R code in the "pandemic response" branch is designed to estimate the structual preference parameters of state govenments. 
Each state faces a trade-off between employment and deaths, and we seek to estimate the relative value of life based on the
observed sequence of lockdown orders. To run this code, first use the file "Covid Data Smoothing". This file runs a smoothing
procedure to recover the infection-spread parameter of a SIRD model. It then estimates a panel relationship between lockdown orders
pandemic spread, and labor market conditions.

Next, use the file CF generator V2 to generate counterfactual paths of COVID-19 spread and employment conditions in each state. The file 
SMM estimator then conductions the minimum-distnace estimation procedure.

