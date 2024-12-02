"""Provide aliases for some common terms.

The modules provides a common list of aliases for different concepts
often used in ESTAMPES.
"""

spectroscopy = {
    # # Vibrational
    # IR
    'Infrared'.casefold(): 'IR',
    'IR'.casefold(): 'IR',
    # VCD
    'VibrationalCircularDichroism'.casefold(): 'VCD',
    'VCD'.casefold(): 'VCD',
    # RS
    'Raman'.casefold(): 'RS',
    'RamanScattering'.casefold(): 'RS',
    'RS'.casefold(): 'RS',
    # ROA
    'RamanOpticalActivity'.casefold(): 'ROA',
    'ROA'.casefold(): 'RS',
    # # Electronic/Resonances
    # OPA
    'OnePhotonAbsorption'.casefold(): 'OPA',
    'Absorption'.casefold(): 'OPA',
    'OPA'.casefold(): 'OPA',
    # OPE
    'OnePhotonEmission'.casefold(): 'OPE',
    'Fluorescence'.casefold(): 'OPE',
    'OPE'.casefold(): 'OPE',
    # OPP
    'OnePhotonPhosphorescence'.casefold(): 'OPP',
    'Phosphorescence'.casefold(): 'OPP',
    'OPP'.casefold(): 'OPP',
    # ECD
    'ElectronicCircularDichroism'.casefold(): 'ECD',
    'CircularDichroism'.casefold(): 'ECD',
    'ECD'.casefold(): 'ECD',
    'CD'.casefold(): 'ECD',
    # CPL
    'CircularPolarizedLuminescence'.casefold(): 'CPL',
    'CPL'.casefold(): 'CPL',
    # CPP
    'CircularPolarizedPhosphorescence'.casefold(): 'CPP',
    'CPP'.casefold(): 'CPP',
    # CPP
    'CircularPolarizedPhosphorescence'.casefold(): 'CPP',
    'CPP'.casefold(): 'CPP',
    # RR
    'ResonanceRaman'.casefold(): 'RR',
    'vRR'.casefold(): 'RR',
    'RR'.casefold(): 'RR',
    # RR
    'ResonanceRamanOptical Activity'.casefold(): 'RROA',
    'vRROA'.casefold(): 'RROA',
    'RROA'.casefold(): 'RROA',
}

level_theory = {
    # # vibrational
    # harmonic
    'Harmonic'.casefold(): 'H',
    'Harm'.casefold(): 'H',
    'Harm.'.casefold(): 'H',
    'H'.casefold(): 'H',
    # anharmonic
    'Anharmonic'.casefold(): 'A',
    'Anharm'.casefold(): 'A',
    'Anharm.'.casefold(): 'A',
    'Anh'.casefold(): 'A',
    'A'.casefold(): 'A',
    # # Electronic
    # pure electronic
    'Electronic'.casefold(): 'E',
    'VerticalExcitation'.casefold(): 'E',
    'PureVE'.casefold(): 'E',
    'E'.casefold(): 'E',
    # vibronic
    'VibrationallyResolvedElectronic'.casefold(): 'VE',
    'Vibronic'.casefold(): 'VE',
    'FranckCondon'.casefold(): 'VE',
    'Franck-Condon'.casefold(): 'VE',
    'FC'.casefold(): 'VE',
    'VE': 'VE'
}

broadening_functions = {
    # Lorentzian broadening
    'lorentzian'.casefold(): 'lorentzian',
    'lorentz'.casefold(): 'lorentzian',
    'lor'.casefold(): 'lorentzian',
    'lor.'.casefold(): 'lorentzian',
    'L'.casefold(): 'lorentzian',
    # Gaussian broadening
    'gaussian'.casefold(): 'gaussian',
    'gauss'.casefold(): 'gaussian',
    'gau'.casefold(): 'gaussian',
    'gau.'.casefold(): 'gaussian',
    'G'.casefold(): 'gaussian',
}
