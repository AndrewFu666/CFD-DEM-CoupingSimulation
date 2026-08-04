#include "les_config.h"

Config::Config(): 
    ReTau(180), 
    Visc(1.0 / 180),
    Fx(1.0),
    dt(1e-3), 
    MaxStp(10000), 
    StatsStp(5000),
    Um(20), A(0.05),
    q(0.5)
{}