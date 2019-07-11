# Tips for hybrid functional calculations
-----------------------

Hybrid functionals, especially HSE06 functional and that with different exchange mixing parameter or screening distance, have
been regularly used for point-defect calculations.

Usually, calculations of hybrid functional are a few hundreds more expensive than the GGA functionals.
Therefore, we need to take a little ingenuity to reduce their computational costs.
For this purpose, we regularly prepare the WAVECAR files obtained using GGA.
(Although we also relax the atomic positions using GGA beforehand in some cases, it would not be good for point-defect
calculations. The reason is because site symmetry of a defect estimated by GGA could be different from that by hybrid functionals.
Furthermore, electronic structures of defects could also be different.)

One can create the INCAR file for generating WAVECAR files using the GGA with the following command, for instance.
```
grep -v LHFCALC INCAR | grep -v ALGO | sed s/"NSW     =  50"/"NSW     =   1"/ > INCAR-pre
```

