#module CellularAutomata
using Distributions
using DelimitedFiles
using Random
using StatsBase

#export colonise,coloniseMultiple,extinction,extinctionNeighbour


struct OccurenceCellularAutomata
    pa::Array
    pa_cart_index::Array
    suitability::Array
    dispersalProbability::Float64
    maxNumberDispersers::Int64
    neighbourSurvivalWeight::Float64
    neighbourDispersalMultiplier::Float64
    meanDispersal::Float64
end

"""
    colonise(ca)

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
function newPos(meanDistance,cartIdx)
    distance = rand(Exponential(meanDistance),1)
    # + 1 ensures dispersal outside of the initial cell
    distance = distance[1] + 0.75
    angle = 360*rand()
    # Remember 1 = y, 2 = x
    x = Int(round(cos(deg2rad(angle))*distance,digits=0)) + cartIdx[2]
    y = Int(round(sin(deg2rad(angle))*distance,digits=0)) + cartIdx[1]
    return(y,x)
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
function coloniseNeighbourWeight(ca::OccurenceCellularAutomata)
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    idxShuffle = sample(collect(1:1:length(dCells)),length(dCells),replace=false)
    for i in idxShuffle
        cellIndex = dCells[i]
        pctColonised = percentColonisedNeighbours(ca.pa,cellIndex,shape)
        # Determine maximum number of dispersers scaled by suitabiility
        maxNumDispersers = ceil(ca.maxNumberDispersers*ca.suitability[i] * (pctColonised * ca.neighbourDispersalMultiplier),digits=0)
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

function percentColonisedNeighbours(pa,idx,shape)
    # Define 3x3 neighbourhood boundaries
    ymin = idx[1]-1
    ymax = idx[1]+1
    xmin = idx[2]-1
    xmax = idx[2]+1
    # Check boundary conditions
    if ymin <= 0
        ymin = 1
    end
    if ymax > shape[2]
        ymax = shape[2]
    end
    if xmin <= 0
        xmin = 1
    end
    if xmax > shape[2]
        xmax = shape[2]
    end
    # Minus 1 to remove the focal cell
    pctColonised =  (sum(pa[ymin:ymax,xmin:xmax]) - 1)/8
    return pctColonised
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
            pctColonised = percentColonisedNeighbours(ca.pa,cellIndex,shp)
            neighboursWeight = pctColonised * ca.neighbourSurvivalWeight
            survivalProbability = ca.suitability[cellIndex] * (1 + neighboursWeight)
            if survivalProbability<1.0
            # Determine survival
                survived = sample([0,1],ProbabilityWeights([1.0-survivalProbability,survivalProbability]))
                if survived === 0
                    ca.pa[cellIndex] = 0
                end
            end
        end
    end
end

function rescale(x,newMin,newMax)
    return (newMax-newMin).*((x.-minimum(x))./(maximum(x)-minimum(x))).+newMin
end

function simulate()
    # Parameters
    landscapeScaleFactors = [0.8]
    dispersalProbability = 0.6
    maxNumberDispersers = [5]
    neighbourSurvivalWeight = [0.4]
    neighbourDispersalMultiplier = [5.0]
    meanDispersalDistance = 1.0
    n_iterations = 100
    for lsf in landscapeScaleFactors
        for mnd in maxNumberDispersers
            for nsw in neighbourSurvivalWeight
                for ndm in neighbourDispersalMultiplier
                    landscape = readdlm("D:/PHDExperimentOutputs/SimLandscapes/suitability/LLM1_suitability.asc",skipstart=6)
                    ls_dimension = size(landscape)
                    pa = ones(ls_dimension)
                    pa_cart_index = CartesianIndices(pa)
                    landscape = rescale(landscape,0,lsf)
                    #testCA = OccurenceCellularAutomata(pa,pa_cart_index,landscape,dispersalProbability,distance,distanceProbabilities,5,0.4,3.0) @ SAc 0.1
                    testCA = OccurenceCellularAutomata(pa,pa_cart_index,landscape,dispersalProbability,mnd,nsw,ndm,meanDispersalDistance)
                    for i in 1:n_iterations
                        #colonise(pa,pa_cart_index,dispersalProbability,distance,distanceProbabilities,ls_dimension)
                        coloniseNeighbourWeight(testCA)
                        #coloniseMultiple(testCA)
                        #extinction(pa,pa_cart_index,landscape)
                        extinctionNeighbour(testCA)
                        # colonise(testCA)
                        # extinction(testCA)
                    end
                    open("D:/test/test_"*string(lsf)*"_"*string(mnd)*"_"*string(nsw)*"_"*string(ndm)*"_"*".asc","w") do io
                        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                        writedlm(io,pa)
                    end
                end
            end
        end
    end
end
#end # Module
@time simulate()
