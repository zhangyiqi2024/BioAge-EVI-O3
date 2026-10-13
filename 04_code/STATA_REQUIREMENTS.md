# Stata requirements

Tested with Stata 18 MP on Windows. Stata itself is proprietary and not bundled.
The analysis code uses standard Stata commands plus estout/esttab, coefplot,
ivreg2 and ranktest. Install these into your licensed Stata environment if absent:

    ssc install estout
    ssc install coefplot
    ssc install ivreg2
    ssc install ranktest

The run log records Stata's header and executed commands. Package versions can
change; the actual resolved ado paths/version headers are captured by the master.
Do not silently update packages in the middle of a reproduction run.
