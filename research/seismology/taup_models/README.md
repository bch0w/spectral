# ak135f_upper_crust

`ak135f_upper_crust.nd`/`.npz` is a custom ObsPy TauP model that, unlike
ObsPy's own `ak135f_no_mud`, it keeps ak135f's
real crustal discontinuities at 10 km and 18 km instead of discarding them:
`ak135f_no_mud` replaces ak135f's entire shallow (<120 km) structure with the
flat, deeper `ak135` crust (Moho at 35 km), whereas this model only strips out
ak135f's thin near-surface water/mud layer (Vp=1.45–1.65 km/s down to 3.3 km)
and replaces it with the upper-crust velocities extended up to the surface,
leaving the 10/18 km crustal structure otherwise intact — so it behaves like
true ak135f without the near-source ray-tracing problems the zero-shear-
velocity mud layer causes. It also places the `mantle` label that TauP uses to
define `moho_branch` (which controls where `Pn`/`Sn` head waves and
`PmP`/`SmS`-type reflections refract) at 18 km, the real velocity
discontinuity, rather than at `ak135f_no_mud`'s 35 km Moho or stock ak135f's
80 km (the base of its 18–80 km crust/mantle gradient zone)
