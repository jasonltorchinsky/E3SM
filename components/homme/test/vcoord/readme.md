# Vertical Coordinate

Stand-alone HOMME utilize a hybrid sigma vertical coordinate, which is of the form
```
p(eta) = A(eta) * p0 + B(eta) * ps
```
where `p` is the pressure at level `eta`, `p0 = 1000 hPa` is a reference pressure, 
`ps` is the pressure at model bottom, and `A` and `B` are coefficients controlling
the rate of transition from terrain-following (`A = 1`, `B = 0`) to pure-pressure
(`A = 0`) vertical coordinate.

# Coefficient Files
These files separate the `A` and `B` coefficients into interface and midpoint files,
with coefficients corresponding the vertical level interfaces and midpoints.

File-pairs are often also labelled with the number of vertical level midpoints.

Different file pairs correspond to different vertical coordinates with different
extents, and often have different naming conventions. Here are my best guesses,
with additional information:
- `sab*`: These are pure-pressure coordinates with a model top nominally at 0 hPa.
- `scream`: Vertical coordinate used in SCREAM. The highest model level midpoint is
at approximately 2.5 hPa.
- `cam`: Vertical coordinate used in CAM. The highest model level midpoint is at
approximately 3.5 hPa.
- `acme` : Vertical coordinate used in EAMxx (?). The highest model level midpoint
is at approximately 0.15 hPa.
- `turbeville` : Vertical coordinate used in numerical experiments for the SCREAM-STRAT project.