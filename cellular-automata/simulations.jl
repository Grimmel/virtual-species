include("ca.jl")

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
runSimulations2()
