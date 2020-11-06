include("ca.jl")
using .CellularAutomata
using DelimitedFiles
using DataFrames

function runSimulations()
    # Parameters
    landscapes = ("LLM","LMS","LLS","LMM","LSS")
    landscapeReplicates = 10
    simulationReplicates = 9 # 9 is the max, 10 total
    speciesList = Dict("ca_species1"=>[0.6,[1],[1]],
                        "ca_species2"=>[0.6,[1,2,3,4],[0.7,0.2,0.07,0.03]],
                        "ca_species3"=>[0.8,[1],[1]],
                        "ca_species4"=>[0.8,[1,2,3,4],[0.7,0.2,0.07,0.03]])

    pa = ones(400,400)
    pa_cart_index = CartesianIndices(pa)
    for species in speciesList
        for landscape in landscapes
            for lsrep in 1:10
                for srep in 0:9
                    landscape = readdlm("D:/PHDExperimentOutputs/SimLandscapes/suitability/"*landscape*string(lsrep)*"_suitability.asc",skipstart=6)
                    landscape = landscape ./100
                    ls_dimension = size(landscape)
                    # Dispersal Parameters
                    dispersalProbability = species[2][1]
                    distance = species[2][2]
                    distanceProbabilities = species[2][3]
                    n_iterations = 200
                    for i in 1:n_iterations
                        colonise(pa,pa_cart_index,dispersalProbability,distance,distanceProbabilities,ls_dimension)
                        extinction(pa,pa_cart_index,landscape)
                    end
                    open("D:/PHDExperimentOutputs/MainSims/"*species[1]*"/Output_maps/abundance/abundance_s1_"*landscape*string(lsrep)*"r"*string(srep)*".tif","w") do io
                        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                        writedlm(io,pa)
                    end
                end
            end
        end
    end
end
function runSimulations2()
    # Parameters
    landscapes = ("LLM","LMS","LLS","LMM","LSS")
    landscapeReplicates = 10
    simulationReplicates = 10
    speciesList = Dict(
                        "ca_species6"=>[0.6,[1,2,3,4],[0.7,0.2,0.07,0.03],1],
                        "ca_species7"=>[0.6,[1],[1],5])
                        #"ca_species8"=>[0.6,[1,2,3,4],[0.7,0.2,0.07,0.03],5])#"ca_species5"=>[0.6,[1],[1],1],
    pa = ones(400,400)
    pa_cart_index = CartesianIndices(pa)
    for species in speciesList
        for landscape in landscapes
            for lsrep in 1:10
                for srep in 0:9
                    suit = readdlm("D:/PHDExperimentOutputs/SimLandscapes/suitability/"*landscape*string(lsrep)*"_suitability.asc",skipstart=6)
                    suit = suit ./100
                    ls_dimension = size(suit)
                    # Dispersal Parameters
                    dispersalProbability = species[2][1]
                    distance = species[2][2]
                    distanceProbabilities = species[2][3]
                    max_dispersers = species[2][4]
                    n_iterations = 200
                    for i in 1:n_iterations
                        if max_dispersers === 1
                            colonise(pa,pa_cart_index,dispersalProbability,distance,distanceProbabilities,ls_dimension)
                        else
                            colonise_multiple(pa,pa_cart_index,suit,dispersalProbability,distance,distanceProbabilities,ls_dimension,max_dispersers)
                        end
                        extinctionNeighbour(pa,pa_cart_index,suit,0.15,ls_dimension)
                    end
                    open("D:/PHDExperimentOutputs/MainSims/"*species[1]*"/Output_maps/abundance/abundance_s1_"*landscape*string(lsrep)*"_r"*string(srep)*".tif","w") do io
                        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                        writedlm(io,pa)
                    end
                end
            end
        end
    end
end
function rescale(x,newMin,newMax)
    return (newMax-newMin).*((x.-minimum(x))./(maximum(x)-minimum(x))).+newMin
end
function rescale(x,newMin,newMax,oldMin,oldMax)
    return (newMax-newMin).*((x.-oldMin)./(oldMax-oldMin)).+newMin
end
function testParams()
    scaleMax = [0.6,0.7,0.8]
    nDisp = [5,2]
    survWeight = [0.5,1.0,1.5]
    meanDisp = [1.0,5.0]
    dispProb = [0.8,0.6]
    for scalem in scaleMax
        suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability238.asc",skipstart=6)
        suit = rescale(suit,0,scalem,0,1)
        for maxDis in nDisp
            for weight in survWeight
                for mdis in meanDisp
                    for dp in dispProb
                        for i in 1:10
                            simName = "sim_"*string(scalem)*"_"*string(maxDis)*"_"*string(weight)*"_"*string(mdis)*"_"*string(dp)*"_rep"*string(i)
                            pa = ones(400,400)
                            paIdx = CartesianIndices(pa)
                            ca = OccurenceCellularAutomata(pa,paIdx,suit,dp,maxDis,weight,1.0,mdis)
                            for j in 1:100
                                coloniseNeighbourWeight(ca)
                                extinctionNeighbour(ca)
                            end
                            open("D:/PHDExperimentOutputs/Transferability/ca_tuning/r1/"*simName*".asc","w") do io
                                write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                                writedlm(io,pa)
                            end
                        end
                    end
                end
            end
        end
    end
end
function testParams2()
    survWeight = [0.2,0.7,1.0,1.5,2.0]
    dispWeight = [0.2,0.4,0.6,0.8,1.0]
    suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability238.asc",skipstart=6)
    suit = rescale(suit,0,0.8,0,1)
    for sweight in survWeight
        for dweight in dispWeight
            for i in 1:10
                simName = "sim_"*string(sweight)*"_"*string(dweight)*"_rep"*string(i)
                pa = ones(400,400)
                paIdx = CartesianIndices(pa)
                ca = OccurenceCellularAutomata(pa,paIdx,suit,0.7,2,sweight,sweight,3.0)
                for j in 1:100
                    coloniseNeighbourWeight(ca)
                    extinctionNeighbour(ca)
                end
                open("D:/PHDExperimentOutputs/Transferability/ca_tuning/r3/"*simName*".asc","w") do io
                    write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                    writedlm(io,pa)
                end
            end
        end
    end
end
function testParams3()
    survWeight = [0.1,0.5,1,3]
    nbR = [1,2,3,4]
    suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability232.asc",skipstart=6)
    suit = rescale(suit,0,0.75,0,1)
    for sweight in survWeight
        for r in nbR
            m = idwDistWeightMatrix(4,1)
            weightMatrix = idwDistWeightMatrix(r,10)
            for i in 1:10
                simName = "ns10sim_"*string(sweight)*"_"*string(r)*"_rep"*string(i)
                pa = ones(400,400)
                paIdx = CartesianIndices(pa)
                ca = OccurenceCellularAutomata(pa,paIdx,suit,2.0,0.4,3,sweight,1.0,weightMatrix,r)
                for j in 1:100
                    coloniseNeighbourWeight(ca)
                    extinctionNeighbourWeight(ca)
                end
                open("D:/PHDExperimentOutputs/Transferability/ca_tuning/r6/"*simName*".asc","w") do io
                    write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                    writedlm(io,pa)
                end
            end
        end
    end
end
function testParams4()
    survWeight = [0.01,0.05]
    nbR = [1,2,3]
    scaleFactor = [0.6,0.75,"na"]
    landscapes = [232,238,218,188]
    for ls in landscapes
        for sf in scaleFactor
            suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability"*string(ls)*".asc",skipstart=6)
            if sf != "na"
                suit = rescale(suit,0,sf,0,1)
            end
            for sweight in survWeight
                for r in nbR
                    weightMatrix = idwDistWeightMatrix(r)
                    for i in 1:10
                        simName = "ls_"*string(ls)*"sim_"*string(sweight)*"_"*string(sf)*"_"*string(r)*"_rep"*string(i)
                        pa = ones(400,400)
                        paIdx = CartesianIndices(pa)
                        ca = OccurenceCellularAutomata(pa,paIdx,suit,2,0.3,4,sweight,1.0,weightMatrix,r)
                        for j in 1:100
                            coloniseNeighbourWeight(ca)
                            extinctionNeighbourWeight(ca)
                        end
                        open("D:/PHDExperimentOutputs/Transferability/ca_tuning/r4/"*simName*".asc","w") do io
                            write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                            writedlm(io,pa)
                        end
                    end
                end
            end
        end
    end
end
function testParamsDisp()
    distance = [1,3,6]
    nDisperser = [1,2,5]
    nbWeight = [0.01,0.5,1.0]
    prob = [0.1,0.5,0.8]
    suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability232.asc",skipstart=6)
    suit = rescale(suit,0,0.75,0,1)
    for d in distance
        for ndisp in nDisperser
            for nb in nbWeight
                for p in prob
                    weightMatrix = idwDistWeightMatrix(1)
                    for i in 1:10
                        simName = "sim_"*string(d)*"_"*string(ndisp)*"_"*string(p)*"_"*string(nb)*"_rep"*string(i)
                        pa = ones(400,400)
                        paIdx = CartesianIndices(pa)
                        ca = OccurenceCellularAutomata(pa,paIdx,suit,d,p,ndisp,0.5,nb,weightMatrix,1)
                        for j in 1:100
                            coloniseNeighbourWeight(ca)
                            extinctionNeighbourWeight(ca)
                        end
                        open("D:/PHDExperimentOutputs/Transferability/ca_tuning/r5/"*simName*".asc","w") do io
                            write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                            writedlm(io,pa)
                        end
                    end
                end
            end
        end
    end
end
function test()
    distance = [1,3,6]
    nDisperser = [1,2,5]
    nbWeight = [0.01,0.5,1.0]
    prob = [0.1,0.5,0.8]
    suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability232.asc",skipstart=6)
    suit = rescale(suit,0,0.75,0,1)
    weightMatrix = idwDistWeightMatrix(1)
    pa = ones(400,400)
    paIdx = CartesianIndices(pa)
    ca = OccurenceCellularAutomata(pa,paIdx,suit,2,0.3,3,0.5,0.1,weightMatrix,1)
    for j in 1:100
        coloniseNeighbourWeight(ca)
        extinctionNeighbourWeight(ca)
    end
    # open("D:/PHDExperimentOutputs/Transferability/ca_tuning/r5/"*simName*".asc","w") do io
    #     write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
    #     writedlm(io,pa)
    # end
end
function testLHS()
    lhsParamSet = readdlm("D:/PHDExperimentOutputs/Transferability/ca_tuning/LHS_1.csv", ',', Float64,skipstart=1)
    for paramSet in lhsParamSet
        print(paramSet)
        param_sample_id = paramSet[1]
        mean_dispersal = paramSet[2]
        dispersal_prob = paramSet[3]
        disp_nb_weight = paramSet[4]
        surv_nb_weight = paramSet[5]
        suit_max = paramSet[6]
        nb_radius = Int(paramSet[7])
        println(nb_radius)
        nb_weight_mat = idwDistWeightMatrix(nb_radius,10)
        max_dispersers = paramSet[8]

        suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability232.asc",skipstart=6)
        suit = rescale(suit,0,suit_max,0,1)
        for i in 1:10
            simName = "sim_"*string(param_sample_id)*"_ls232_rep"*string(i)
            pa = ones(400,400)
            paIdx = CartesianIndices(pa)
            ca = OccurenceCellularAutomata(pa,paIdx,suit,mean_dispersal,dispersal_prob,max_dispersers,surv_nb_weight,disp_nb_weight,nb_weight_mat,nb_radius)
            for j in 1:100
                coloniseNeighbourWeight(ca)
                extinctionNeighbourWeight(ca)
            end
            open("D:/PHDExperimentOutputs/Transferability/ca_tuning/lhs1/"*simName*".asc","w") do io
                write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                writedlm(io,pa)
            end
        end
    end
end

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
testLHS()
