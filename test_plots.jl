
# @benchmark Posc_vacuum_planewave( mu, mu, natural(30km)/200MeV )

# Posc_vacuum_planewave( mu, mu, natural(30km)/200MeV )
# Posc_vacuum_planewave( mu, ele, natural(30km)/200MeV )
# Posc_vacuum_planewave( mu, tau, natural(30km)/200MeV )


# check solar oscillations (@ KamLAND )
begin 
    α = ele
    L = (50:0.2:200) * u"km"
    E = 5MeV

    lw = 1.5
    β = ele
        P = Posc_vacuum_planewave.( α, β, natural.(L) ./ E )
        plot( L, P; lw, xlabel="baseline L", label="P($α → $β)" )
    #
    
    β = mu
        P .= Posc_vacuum_planewave.( α, β, natural.(L) ./ E )
        plot!( L, P; lw, xlabel="baseline L", label="P($α → $β)" )
    #

    β = tau
        P .= Posc_vacuum_planewave.( α, β, natural.(L) ./ E )
        plot!( L, P; lw, xlabel="baseline L", label="P($α → $β)" )
    #

    plot!( title="\nE = $E", ylims=(-0.05,1.05), size=(500, 250) )

end


# check atmospheric oscillations 
begin 

    α = mu

    L = 30km
    E = exp10.(-1.5:0.001:3) * GeV

    lw = 1.5
    β = ele
        P = Posc_vacuum_standard.( α, β, natural.(L) ./ E )
        plot( E, P; lw, xlabel="energy E", label="P($α → $β)" )
    #
    
    β = mu
        P .= Posc_vacuum_standard.( α, β, natural.(L) ./ E )
        plot!( E, P; lw, label="P($α → $β)" )
    #

    β = tau
        P .= Posc_vacuum_standard.( α, β, natural.(L) ./ E )
        plot!( E, P; lw, label="P($α → $β)" )
    #

    xt = [10MeV, 30MeV, 100MeV, 300MeV, 1GeV, 3GeV, 10GeV]
    sxt = ["10MeV", "30MeV", "100MeV", "300MeV", "1GeV", "3GeV", "10GeV"]

    plot!(
        xscale=:log10, xticks=( xt, sxt ),
        title="\nL = $L",
        size=(500, 250)
    )
    
end