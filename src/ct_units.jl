"""

Physical units in SI units

"""

# SI units
const K = 1.0         # Kelvin
const J = 1.0         # Joule
const A = 1.0         # Ampere
const V = 1.0         # Volt
const m = 1.0         # meter
const s = 1.0         # seconds
const C = A * s       # Coulomb
const kg = 1.0         # kilogram
const Hz = 1 / s       # Hertz
const kHz = 1.0e3 * Hz # Kilohertz
const W = 1.0        # Watt
const kW = 1.0e3 * W  # KiloWatt


# scaled SI units
const cm = 1.0e-2 * m
const mm = 1.0e-3 * m
const μm = 1.0e-6 * m
const nm = 1.0e-9 * m

const ms = 1.0e-3 * s
const μs = 1.0e-6 * s
const ns = 1.0e-9 * s
const ps = 1.0e-12 * s

# not constant: depends on the given value of q
global eV = q * V
