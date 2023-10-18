# TurbulenceModels

Custom OpenFOAM v2306 TurbulenceModels library version.
Differences between the vanilla and the present version:
* S3PQR family [1] models included within LES.
* S3PQR family [1] models included within DDES.

** DDES solver modified in order to include advanced turbulence models.
* Delta_lsq [2] included within LESDeltas
* Deck (2020) [3] shielding function included within DDES.

[1] Trias, F. X., Folch, D., Gorobets, A., & Oliva, A. (2015). Building proper invariants for eddy-viscosity subgrid-scale models. Physics of Fluids, 27(6). https://doi.org/10.1063/1.4921817

[2] Trias, F. X., Gorobets, A., Silvis, M. H., Verstappen, R. W. C. P., & Oliva, A. (2017). A new subgrid characteristic length for turbulence simulations on anisotropic grids. Physics of Fluids, 29(11). https://doi.org/10.1063/1.5012546

[3] Deck, S., & Renard, N. (2020). Towards an enhanced protection of attached boundary layers in hybrid RANS/LES methods. Journal of Computational Physics, 400, 108970. https://doi.org/10.1016/j.jcp.2019.108970
