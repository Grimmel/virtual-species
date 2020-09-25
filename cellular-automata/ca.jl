using Distributions
using DelimitedFiles
using Random
using StatsBase

function colonise(pa,pa_cart_index,dispersal_prob,distances,dist_prob,shape)
    for i in pa_cart_index
        if pa[i] === 1.0
            disperse = rand(Bernoulli(dispersal_prob),1)
            if disperse[1] == true
                distance = sample(distances,ProbabilityWeights(dist_prob))
                angle = 360*rand()
                # 1 = y, 2 = x
                x = Int(round(cos(deg2rad(angle))*distance,digits=0)) + i[2]
                y = Int(round(sin(deg2rad(angle))*distance,digits=0)) + i[1]
                if x>=1 && x <= shape[2] && y >=1 && y<=shape[1]
                    pa[y,x] = 1.0
                end
            end
        end
    end
end

function colonise_multiple(pa,pa_cart_index,landscape,dispersal_prob,distances,dist_prob,shape,max_dispersers)
    for i in pa_cart_index
        if pa[i] === 1.0
            # Determine maximum number of dispersers scaled by suitabiility
            maxNumDispersers = ceil(max_dispersers*landscape[i],digits=0)
            potentialDispersers = collect(0:1:maxNumDispersers)
            # Weights are the inveverse of number of dispersers,
            # normalised to sum to 1
            sampleWeights = potentialDispersers.^-1
            replace!(sampleWeights,Inf=>0)
            # Ensure 0 number of dispersers are weighted equal to 1-dispersal_prob
            sampleWeights[1] = sum(sampleWeights/dispersal_prob) * (1-dispersal_prob)
            sampleWeights = sampleWeights./sum(sampleWeights)
            numberDispersers = sample(potentialDispersers,ProbabilityWeights(sampleWeights))
            for j in 1:numberDispersers
                distance = sample(distances,ProbabilityWeights(dist_prob))
                angle = 360*rand()
                # Remember 1 = y, 2 = x
                x = Int(round(cos(deg2rad(angle))*distance,digits=0)) + i[2]
                y = Int(round(sin(deg2rad(angle))*distance,digits=0)) + i[1]
                if x>=1 && x <= shape[2] && y >=1 && y<=shape[1]
                    pa[y,x] = 1.0
                end
            end
        end
    end
end

function extinction(pa,pa_cart_index,landscape)
    for i in pa_cart_index
        if pa[i] === 1.0
            suitabilty = landscape[i]
            survived = rand(Bernoulli(suitabilty),1)
            if survived[1] == false
                pa[i] = 0
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

function extinctionNeighbour(pa,pa_cart_index,landscape,maxWeight,shp)
    # Randomise cell index to prevent bias from the neighbourhood weighting
    # Alternative could be to copy the the pa array.
    nCells = shp[1]*shp[2]
    idxShuffle = sample(collect(1:1:nCells),nCells,replace=false)
    for i in idxShuffle
        cellIndex = pa_cart_index[i]
        if pa[i] === 1.0
            # Determine the neighbourhood weighted survival probability
            pctColonised = percentColonisedNeighbours(pa,cellIndex,shp)
            neighboursWeight = pctColonised * maxWeight
            survivalProbability = landscape[cellIndex] * (1 + neighboursWeight)
            if survivalProbability>1
                survivalProbability = 1
            end
            # Determine survival
            survived = rand(Bernoulli(survivalProbability),1)
            if survived[1] == false
                pa[cellIndex] = 0
            end
        end
    end
end

function simulate()
    landscape = readdlm("D:/PHDExperimentOutputs/SimLandscapes/suitability/LLM1_suitability.asc",skipstart=6)
    landscape = landscape ./100
    ls_dimension = size(landscape)
    pa = ones(ls_dimension)
    pa_cart_index = CartesianIndices(pa)
    # Dispersal Parameters
    dispersalProbability = 0.6
    distance = [1,2,3]
    distanceProbabilities = [0.85,0.14,0.1]
    n_iterations = 200
    for i in 1:n_iterations
        colonise_multiple(pa,pa_cart_index,landscape,dispersalProbability,distance,distanceProbabilities,ls_dimension,5)
        #extinction(pa,pa_cart_index,landscape)
        extinctionNeighbour(pa,pa_cart_index,landscape,0.5,ls_dimension)
    end
    open("D:/test/julia_test.asc","w") do io
        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
        writedlm(io,pa)
    end
end
#@time simulate()
