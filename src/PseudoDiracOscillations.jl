module PseudoDiracOscillations

# using StaticArrays
using Rotations

using NaturallyUnitful
using UnitfulAstro
include( "units.jl" )


@enum Flavor ele=1 mu=2 tau=3
Base.to_index( α::Flavor ) = Base.to_index( Int(α) )

function Posc_vacuum_standard( α, β, LoE )

    # NuFit w/o SK; NO 

        Δmsq_21 = 1e-5eV^2 * 7.41 # ± (0.21, 0.20)
        Δmsq_3l = 1e-3eV^2 * 2.511 # ± (0.027, 0.027)

        Δmsq_31 = Δmsq_3l
        Δmsq_32 = Δmsq_31 - Δmsq_21 

        θ12 = 33.66° # ± ( 0.73, 0.70 )
        θ23 = 49.1° # ± ( 1.0, 1.3 )
        θ13 = 8.54° # ± ( 0.11, 0.11 )
        δCP = 197° # ± (41, 25)
    
    # 

    # PMNS matrix
        R23 = RotX( θ23 )
        R13 = RotMatrix{3}(
            [
                cos(θ13) 0 sin(θ13) * exp( -im * δCP ); 
                0 1 0;
                -sin(θ13) * exp( im * δCP ) 0 cos(θ13);
            ]
        )
        R12 = RotZ( θ12 )
        U = (R23 * R13 * R12)
    # 

    z(i, j) = U[α,i] * conj(U[β,i]) * conj(U[α,j]) * U[β,j]
    
    # i < j
    ij_iter = ( (1,3), (2,3), (1,2) )
    Δmsq_ij_iter = ( Δmsq_31, Δmsq_32, Δmsq_21 )

    P = float( α==β )
    for (ij, Δmsq_ij) in zip( ij_iter, Δmsq_ij_iter )

        re_ij, im_ij = reim( z( ij... ) )

        P += -4 * re_ij * sin( -Δmsq_ij * LoE / 4 )^2
        P +=  2 * im_ij * sin( -Δmsq_ij * LoE / 2 )

    end
    return P 
end

function Posc_vacuum_PD( α, β, LoE, δmsq::AbstractArray )

    # NuFit w/o SK; NO 

        θ12 = 33.66° # ± ( 0.73, 0.70 )
        θ23 = 49.1° # ± ( 1.0, 1.3 )
        θ13 = 8.54° # ± ( 0.11, 0.11 )
        δCP = 197° # ± (41, 25)
    
    # 

    # PMNS matrix
        R23 = RotX( θ23 )
        R13 = RotMatrix{3}(
            [
                cos(θ13) 0 sin(θ13) * exp( -im * δCP ); 
                0 1 0;
                -sin(θ13) * exp( im * δCP ) 0 cos(θ13);
            ]
        )
        R12 = RotZ( θ12 )
        U = (R23 * R13 * R12)
    # 

    # assuming standard oscillations average out
    return sum(
        j -> abs2( U[α,j] * U[β,j] ) * cos( δmsq[j] * LoE / 4 )^2,
        1:3
    ) 
    
end

export Flavor, mu, ele, tau
export Posc_vacuum_standard
export Posc_vacuum_PD

# export units
units_to_export = (
    :c,
    :m, :s, :eV,
    :rad, :°, :sr,
    :c0, #:q, :ħ,
    :pc, #:AU, :ly, 
    :MeV, :GeV, :TeV, :PeV, :EeV,
    :kpc, :Gpc, :cm, :km,
    :mb,
    :Gauss, :μG,
)
for x in units_to_export
    @eval export $x
end

export @u_str, natural

end # module PseudoDiracOscillations