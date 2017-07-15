class Procedure:
    NPT = 'npt'
    NVT_CV = 'npt-cv'
    NVT_VACUUM = 'nvt-vacuum'
    NVT_SLAB = 'nvt-slab'
    NPT_BINARY_SLAB = 'npt-binary-slab'
    choices = [NPT, NVT_CV, NVT_VACUUM, NVT_SLAB, NPT_BINARY_SLAB]
    T_RELEVANT = choices
    P_RELEVANT = [NPT, NVT_CV, NPT_BINARY_SLAB]

    prior = {
        NVT_CV : NPT,
    }
