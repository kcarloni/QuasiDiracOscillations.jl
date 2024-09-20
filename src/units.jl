
import NaturallyUnitful: c, natural, unnatural, uconvert
# using UnitfulAstro
import Unitful: @u_str, @unit, unit, Quantity
import Unitful: rad, °, sr

unitful_to_import = (
    :m, :s, :eV,
    :c0, :q, :ħ,
    :yr,
    :b,
    :Gauss
)
for u in unitful_to_import
    @eval import Unitful: $u
end

import UnitfulAstro: AU, ly, pc #, V_mag 
# const c_light = c0

const cm = u"cm"
const km = u"km"

const MeV = u"MeV"
const GeV = u"GeV"
const TeV = u"TeV"
const PeV = u"PeV"
const EeV = u"EeV"

const kpc = u"kpc"
const Mpc = u"Mpc"
const Gpc = u"Gpc"

const mb = u"mb"

# const μG = u"μGauss"


Base.:*( 
    num::Quantity{T1, D, U}, 
    denom::Quantity{T2, D, U}) where {T1, T2, U, D} = 
        num.val * denom.val * unit(num) * unit(denom)
    
Base.:*( 
    num::Quantity{T1, D, U1}, 
    denom::Quantity{T2, D, U2}) where {T1, T2, U1, U2, D} = 
        uconvert(unit(denom), num) * denom

Base.:/( num::Quantity{T1, D, U}, denom::Quantity{T2, D, U}) where {T1, T2, U, D} = num.val / denom.val
Base.:/( num::Quantity{T1, D, U1}, denom::Quantity{T2, D, U2}) where {T1, T2, U1, U2, D} = uconvert(unit(denom), num) / denom
