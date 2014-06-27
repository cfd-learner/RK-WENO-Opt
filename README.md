# Script 

(just change p and s)

    x = load('Eigenvalues_CRWENO5.txt');
    z = x(:,1) + 1i*x(:,2);
    p = 4;
    s = 5;
    [h, poly_coeff] = opt_poly_bisect(z,s,p,'chebyshev','do_plot',true)
    rk = rk_opt(s,p,'erk','ssp','poly_coeff_ind',((p+1):s),'poly_coeff_val',poly_coeff((p+2):s+1),'display','iter')
    plot(real(h*z),imag(h*z),'o')
    hold on
    plotstabreg_func(poly_coeff,[1])

# Results

Values of h/s

## CRWENO5

### p = 4

s  | h/s   
--- | ------
4  | 0.2585
5  | 0.2586
6  | 0.2758
7  | 0.2888
8  | 0.2983
9  | 0.3064
10 | 0.3130
