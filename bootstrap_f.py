def bootstrap_f(tb18v, tb37v):
    """Comiso ice concentration algorithm
    """
    """bootstrap winter tie points (Comiso, 1997)"""
    tiepts = [202.00, 182.00, 246.00, 130.00, 170.00, 234.00, 179.00,\
                  218.00, 253.00, 0.0, 0.0, 0.0]

    tw18v = tiepts[7 - 1]
    tw37v = tiepts[1 - 1]
    tfy18v = tiepts[9 - 1]
    tfy37v = tiepts[3 - 1]
    tmy18v = tiepts[8 - 1]
    tmy37v = tiepts[2 - 1]


    af = (tfy37v - tmy37v)/(tfy18v - tmy18v)
    bf = (tmy37v - af*tmy18v)
    if (tb18v == tw18v): 
        tw18v=tw18v-0.05
    qf = (tb37v - tw37v)/(tb18v - tw18v)
    wf = (tw37v - qf*tw18v)
    ti18vf = (bf - wf)/(qf - af)
    cf = (tb18v - tw18v)/(ti18vf - tw18v)
    cf=min(cf,1.0)
    cf=max(cf,0.0)
    return cf

