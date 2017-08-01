class Procedure:
    NPT = 'npt'
    NVT_CV = 'nvt-cv'
    NVT_VISCOSITY = 'nvt-viscosity'
    NVT_VACUUM = 'nvt-vacuum'
    NVT_SLAB = 'nvt-slab'
    NPT_BINARY_SLAB = 'npt-binary-slab'
    choices = [NPT, NVT_CV, NVT_VISCOSITY, NVT_VACUUM, NVT_SLAB, NPT_BINARY_SLAB]
    T_RELEVANT = choices
    P_RELEVANT = [NPT, NVT_CV, NVT_VISCOSITY, NPT_BINARY_SLAB]

    prior = {
        NVT_CV : NPT,
        NVT_VISCOSITY: NPT,
    }
