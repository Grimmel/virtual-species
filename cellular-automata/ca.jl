module CellularAutomata
using Distributions
using Random
using StatsBase

export colonise,coloniseMultiple,coloniseNeighbourWeight,extinction,extinctionNeighbour,OccurenceCellularAutomata


struct OccurenceCellularAutomata
    pa::Matrix{Float64}
    pa_cart_index::Matrix{CartesianIndex{2}}
    suitability::Matrix{Float64}
    dispersalProbability::Float64
    maxNumberDispersers::Int64
    neighbourSurvivalWeight::Float64
    neighbourDispersalMultiplier::Float64
    meanDispersal::Float64
end


function newPos(meanDistance,cartIdx)
    distance = rand(Exponential(meanDistance),1)
    # + 1 ensures dispersal outside of the initial cell
    distance = distance[1] + 0.75
    angle = 360.0*rand()
    # Remember 1 = y, 2 = x
    x = Int(round(cos(deg2rad(angle))*distance,digits=0)) + cartIdx[2]
    y = Int(round(sin(deg2rad(angle))*distance,digits=0)) + cartIdx[1]
    coord = (y,x)
    return(coord)
end
function selectProportion(pa,caIndex,dispersalProbability)
    nPresences = Int(sum(pa))
    total = Array{CartesianIndex}(undef, nPresences)
    counter = 1
    for idx in caIndex
        if pa[idx] === 1.0
            total[counter] = idx
            counter+=1
        end

    end
    numberDispersing = sample(total,rand(Binomial(nPresences,dispersalProbability)))
end
"""
    colonise(ca)

Select a single cell to colonise

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
function colonise(ca::OccurenceCellularAutomata)
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    for i in dCells
        newXY = newPos(ca.meanDispersal,i)
        if newXY[2]>=1 && newXY[2] <= shape[2] && newXY[1] >=1 && newXY[1]<=shape[1]
            ca.pa[newXY[1],newXY[2]] = 1.0
        end
    end
end

"""
    coloniseSuitWeight(ca)

Select a number of dispersing cells and randomly selects how many cells each can
colonise. A maximum number of dispersers is needed to be defined, with the
actual number of dispersers weighted by suitability.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
function coloniseSuitWeight(ca::OccurenceCellularAutomata)
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    for i in dCells
        # Determine maximum number of dispersers scaled by suitabiility
        maxNumDispersers = ceil(ca.maxNumberDispersers*ca.suitability[i],digits=0)
        potentialDispersers = collect(1:1:maxNumDispersers)
        # Weights are the inveverse of number of dispersers,
        # normalised to sum to 1
        sampleWeights = potentialDispersers.^-1
        replace!(sampleWeights,Inf=>0)
        sampleWeights = sampleWeights./sum(sampleWeights)
        numberDispersers = sample(potentialDispersers,ProbabilityWeights(sampleWeights))
        for j in 1:numberDispersers
            newXY = newPos(ca.meanDispersal,cellIndex)
            if newXY[2]>=1 && newXY[2] <= shape[2] && newXY[1] >=1 && newXY[1]<=shape[1]
                ca.pa[newXY[1],newXY[2]] = 1.0
            end
        end
    end
end
"""
    coloniseNeighbourWeight(ca)

For each occupied cell in the CA, determine the if that cell will disperse and
use an exponential kernel to determine where it will colonise. Distance is
sampled randomly from an exponential distribution, while direction is sampled
from a uniform distribution between 0-360.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
function coloniseNeighbourWeight(ca::OccurenceCellularAutomata)
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    idxShuffle = sample(collect(1:1:length(dCells)),length(dCells),replace=false)
    for i in idxShuffle
        cellIndex = dCells[i]
        neighbourWeight = neighbourHoodWeight(ca.pa,ca.suitability,cellIndex,shape)
        # Determine maximum number of dispersers scaled by suitabiility
        maxNumDispersers = ceil(ca.maxNumberDispersers*ca.suitability[i] * (1+(neighbourWeight * ca.neighbourDispersalMultiplier)),digits=0)
        potentialDispersers = collect(0:1:maxNumDispersers)
        # Weights are the inveverse of number of dispersers,
        # normalised to sum to 1
        sampleWeights = potentialDispersers.^-1
        replace!(sampleWeights,Inf=>0)
        # Ensure 0 number of dispersers are weighted equal to 1-dispersal_prob
        sampleWeights[1] = sum(sampleWeights/ca.dispersalProbability) * (1-ca.dispersalProbability)
        sampleWeights = sampleWeights./sum(sampleWeights)
        numberDispersers = sample(potentialDispersers,ProbabilityWeights(sampleWeights))
        for j in 1:numberDispersers
            newXY = newPos(ca.meanDispersal,cellIndex)
            if newXY[2]>=1 && newXY[2] <= shape[2] && newXY[1] >=1 && newXY[1]<=shape[1]
                ca.pa[newXY[1],newXY[2]] = 1.0
            end
        end
    end
end

function neighbourHoodWeight(pa,suit,idx,shape)
    # Define 3x3 neighbourhood boundaries
    ymin = idx[1]-1
    ymax = idx[1]+1
    xmin = idx[2]-1
    xmax = idx[2]+1
    # Check boundary conditions
    if ymin <= 0
        ymin = 1
    end
    if ymax > shape[1]
        ymax = shape[1]
    end
    if xmin <= 0
        xmin = 1
    end
    if xmax > shape[2]
        xmax = shape[2]
    end
    oc = pa[ymin:ymax,xmin:xmax]
    su = suit[ymin:ymax,xmin:xmax]
    pctColonised = 0
    for i in eachindex(oc)
        pctColonised = pctColonised+(oc[i]*su[i])
    end
    pctColonised = (pctColonised-(pa[idx]*suit[idx]))/8
    return pctColonised
end

function extinction(ca::OccurenceCellularAutomata)
    for idx in ca.pa_cart_index
        if ca.pa[idx] === 1.0
            survived = rand(Bernoulli(ca.suitability[idx]),1)
            if survived[1] == false
                ca.pa[idx] = 0
            end
        end
    end
end
function extinctionNeighbour(ca::OccurenceCellularAutomata)
    # Randomise cell index to prevent bias from the neighbourhood weighting
    # Alternative could be to copy the the pa array.
    shp = size(ca.pa)
    nCells = shp[1]*shp[2]
    idxShuffle = sample(collect(1:1:nCells),nCells,replace=false)
    for i in idxShuffle
        cellIndex = ca.pa_cart_index[i]
        if ca.pa[i] === 1.0
            # Determine the neighbourhood weighted survival probability
            neighbourWeight = neighbourHoodWeight(ca.pa,ca.suitability,cellIndex,shp)
            survWeight = neighbourWeight * ca.neighbourSurvivalWeight
            # Logistic function to scale survival proba
            sf = 1/(1+MathConstants.e^(-10*(ca.suitability[cellIndex]-(1-survWeight))))
            survivalProbability = ca.suitability[cellIndex] + (1-ca.suitability[cellIndex])*sf
            # Determine survival
            survived = sample([0,1],ProbabilityWeights([1.0-survivalProbability,survivalProbability]))
            if survived === 0
                ca.pa[cellIndex] = 0
            end
        end
    end
end

function rescale(x,newMin,newMax)
    return (newMax-newMin).*((x.-minimum(x))./(maximum(x)-minimum(x))).+newMin
end

end # Module
