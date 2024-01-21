lower = c(1)
upper = c(100)

theta = c(0.1,0.1,0.1)
S = BDI.S

Q = get_coffin_matrix(lower,upper,theta,S,"BDI")
Qs = as(Q,"sparseMatrix")

x = index_to_state(c(1),lower,upper)
x
BDI_get_rates(x,theta)
Q

v = matrix(rep(0,nrow(Q)),nrow=1)
obs_index = 4
v[obs_index] = 1

microbenchmark(
    {expm(Q)},
    {SS_exp_Q(Q,prec=1e-15)},
    {SS_exp_Q(Qs,prec=1e-15)},
    {v_exp_Q(v,Q,prec=1e-15)},
    {v_exp_Q(v,Qs,prec=1e-15)},
    unit='micro'
)


