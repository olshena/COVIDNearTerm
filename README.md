# COVIDNearTerm
CovidNearTerm is a bootstrap-based method based on an autoregressive model to estimate at the county level the expected number of COVID-19 patients that will hospitalized 2-4 weeks into the future.  It is based on the work of researchers at UCSF (Adam Olshen), Stanford (Kristopher Kapphahn, Ariadna Garcia, Isabel Wang and Manisha Desai) and Memorial Sloan-Kettering (Mithat Gonen).

Our projections are based on the call simulateAR(vec,wsize=14,method="unweighted",pdays=28,nsim=10000,seed=12345,output_type="predictions",lambda=seq(0,1,0.05))
