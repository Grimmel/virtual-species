include("ca.jl")
using .CellularAutomata
using DelimitedFiles
using DataFrames
using CSV
using GaussianRandomFields
using FileIO
using Colors,Images,FixedPointNumbers
function idwDistWeightMatrix(r,scale)
    w = 2*r+1
    m = ones(w,w)
    mIdx = CartesianIndices(m)
    for idx in mIdx
        #print(idx[1])
        m[idx] = sqrt(((idx[1]-(r+1))^2)+((idx[2]-(r+1))^2))
    end
    m = (m.^-1)/scale
    m[r+1,r+1] = 0.0
    return(m)
end

function rescale(x,newMin,newMax,oldMin,oldMax)
    return (newMax-newMin).*((x.-oldMin)./(oldMax-oldMin)).+newMin
end
struct NegativeInteraction
    species1
    species2
    nicheOverlap
    competition # 1-competition for species 2
end
# struct interactions
#     species
# end
function createInteraction(species)
end
#
# inter1 = Interaction('ca1','ca2',0.5,0.3)
# inter2 = Interaction('ca1','ca3',0.2,0.7)
# inter3 = Interaction('ca2','ca3',0.3,0.5)
# pa1 * pa2 gives highest values where HS of 1 overlap
# Multiplying by nicheOverlap reduces the effect size depending on how much
# each species uses the same resources
# Multiplying by competiton scales competition effect depending on the
# competetive edge of that species
# Multiply by occurence of SP 2 to remove effect of
function simulate()
    suit1 = readdlm("D:/species1.asc",skipstart=6)
    suit2 = readdlm("D:/species2.asc",skipstart=6)
    pa1 = ones(400,400)
    paIdx1 = CartesianIndices(pa1)
    pa2 = ones(400,400)
    paIdx2 = CartesianIndices(pa2)
    weight = idwDistWeightMatrix(2,10)
    ca1 = OccurrenceCellularAutomataNB(deepcopy(pa1),paIdx1,deepcopy(suit1),1,0.8,3,1.0,1.0,weight,2)
    ca2 = OccurrenceCellularAutomataNB(deepcopy(pa2),paIdx2,deepcopy(suit2),1,0.8,3,1.0,1.0,weight,2)
    nicheOverlap = 0.8
    competition = 0.7
    y = zeros(Float64,400,400,200)
    for i in 1:200
        compEffect1_2 = (suit2 .* nicheOverlap.*competition .* ca2.pa).+1
        compEffect2_1 = (suit1 .* nicheOverlap.*(1-competition) .* ca1.pa).+1
        # ca1.suitability = suit1 ./ compEffect1_2
        # ca2.suitability = suit2 ./ compEffect2_1
        coloniseNeighbourWeight(ca1)
        extinctionNeighbourWeight(ca1)
        coloniseNeighbourWeight(ca2)
        extinctionNeighbourWeight(ca2)
    end
    open("D:/sp1_int2.asc","w") do io
        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
        writedlm(io,ca1.pa)
    end
    open("D:/sp2_int2.asc","w") do io
        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
        writedlm(io,ca2.pa)
    end
end
# function interaction(sp1::OccurrenceCellularAutomata,sp2::OccurrenceCellularAutomataNB,overlap::Float64,edge::Float64)
#     e12 = (overlap * edge) .* sp1.suit1
#     e21 = (overlap * (1-edge)) .* sp1.suit1
#     return(e12,e21)
# end
function simulate2()
    suit1 = readdlm("D:/species1.asc",skipstart=6)
    suit2 = readdlm("D:/species2.asc",skipstart=6)
    suit3 = readdlm("D:/species3.asc",skipstart=6)
    pa1 = ones(400,400)
    paIdx1 = CartesianIndices(pa1)
    pa2 = ones(400,400)
    paIdx2 = CartesianIndices(pa2)
    pa3 = ones(400,400)
    paIdx3 = CartesianIndices(pa3)
    weight = idwDistWeightMatrix(2,10)
    ca1 = OccurrenceCellularAutomataNB(deepcopy(pa1),paIdx1,deepcopy(suit1),1,0.8,2,1.0,1.0,weight,2)
    ca2 = OccurrenceCellularAutomataNB(deepcopy(pa2),paIdx2,deepcopy(suit2),1,0.8,2,1.0,1.0,weight,2)
    ca3 = OccurrenceCellularAutomataNB(deepcopy(pa3),paIdx3,deepcopy(suit3),1,0.8,2,1.0,1.0,weight,2)
    s12 = Dict("Species"=>"Species1","Competitor"=>"Species2","Overlap"=>0.5,"Edge"=>0.3)
    s23 = Dict("Species"=>"Species2","Competitor"=>"Species3","Overlap"=>0.9,"Edge"=>0.3)
    richness = zeros(Float64,400,400,200)
    for i in 1:200
        e12 = (s12["Overlap"] * s12["Edge"]) .* suit1
        e21 = (s12["Overlap"] * (1 - s12["Edge"])) .* suit2
        e23 = (s23["Overlap"] * s23["Edge"]) .* suit2
        e32 = (s23["Overlap"] * (1 - s23["Edge"])) .* suit3
        ca1.suitability = suit1.*(1 .- e21 .* ca2.pa)
        ca3.suitability = suit3 .* (1 .- e23 .* ca2.pa)
        ca2.suitability = suit2 .* (1 .- ((e12 .* ca1.pa).+ (e32 .* ca3.pa)))

        coloniseNeighbourWeight(ca1)
        extinctionNeighbourWeight(ca1)
        coloniseNeighbourWeight(ca2)
        extinctionNeighbourWeight(ca2)
        coloniseNeighbourWeight(ca3)
        extinctionNeighbourWeight(ca3)
        richness[:,:,i] = (ca1.pa .+ ca2.pa .+ ca3.pa)./3
    end
    open("D:/sp1_int2.asc","w") do io
        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
        writedlm(io,ca1.pa)
    end
    open("D:/sp2_int2.asc","w") do io
        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
        writedlm(io,ca2.pa)
    end
    open("D:/sp3_int3.asc","w") do io
        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
        writedlm(io,ca3.pa)
    end
    save("D:/rich.gif",richness,fps=10)
end
simulate2()
