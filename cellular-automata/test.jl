include("ca.jl")
using .CellularAutomata
using Plots,Distributions
function test_neighbourhoodWeight()
    pa = [1.0 1.0 0.0;0.0 1.0 0.0; 1.0 0.0 0.0]
    suit= [0.8 0.4 0.4;0.1 0.8 0.7;0.9 0.7 0.3]
    paidx = CartesianIndices(pa)
    expected22 = (suit[paidx[1]]+suit[paidx[4]]+suit[paidx[3]])/8
    expected11 = (suit[4]+suit[5])/8
    expected23 = (suit[4]+suit[5])/8
    wt22 = CellularAutomata.neighbourHoodWeight(pa,suit,paidx[5],size(pa))
    wt11 = CellularAutomata.neighbourHoodWeight(pa,suit,paidx[1],size(pa))
    wt23 = CellularAutomata.neighbourHoodWeight(pa,suit,paidx[8],size(pa))
    return((wt22-expected22),(wt11-expected11),(wt23-expected23))
end
function test_newPos(dispersalDistance,y,x)
    p = zeros(21,21)
    initx = y
    inity = x
    for i in 1:1000000
        pos = CellularAutomata.newPos(dispersalDistance,CartesianIndex(inity,initx))
        if pos[1]>=1 && pos[1]<=21 && pos[2]>=1 && pos[2]<=21
            p[pos[2],pos[1]] +=1
        end
    end
    return(p)
end
function test_multiNewPos()
    t1 = test_newPos(5,11,11)
    t2 = test_newPos(5,5,16)
    t3 = test_newPos(2,11,11)
    t4 = test_newPos(2,5,11)
    plot(heatmap(t1,c = :viridis),heatmap(t2,c = :viridis),heatmap(t3,c = :viridis),heatmap(t4,c = :viridis))
end
function test_coloniseNBW()
    pa = [1.0 1.0 0.0;1.0 1.0 0.0; 1.0 1.0 0.0]
    suit= [0.8 0.4 0.4;0.1 0.8 0.7;0.9 0.7 0.3]
    paidx = CartesianIndices(pa)
    neighbourDispersalWeight=0.5
    maxdis = 5.0
    neighbourWeight = CellularAutomata.neighbourHoodWeight(pa,suit,paidx[5],size(pa))
    dispWeight = neighbourWeight * neighbourDispersalWeight
    # Determine maximum number of dispersers scaled by suitabiility
    dispersalMultiplier = 1/(1+MathConstants.e^(-10*(suit[paidx[5]]-(1-dispWeight))))
    maxNumDispersers = (suit[paidx[5]] + (1-suit[paidx[5]])*dispersalMultiplier)*maxdis
    numberDispersers = rand(Poisson(maxNumDispersers))
    println(maxNumDispersers)
    println(neighbourWeight)
    println(dispWeight)
end
function test_MultiColNBW()

end
