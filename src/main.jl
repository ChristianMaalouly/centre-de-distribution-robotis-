include("lecture.jl")
include("problematique_1.jl")
include("dantzig-wolfe.jl")

data_files = ["Data_test_N10_R10_O10_RS7.txt", "Data_test_N12_R12_O12_RS8.txt", "Data_test_N5_R2_O3_RS2.txt", "Data_test_N5_R3_O3_RS5.txt",
    "Data_test_N5_R4_O3_RS2.txt", "Data_test_N7_R4_O6_RS7.txt", "Data_test_N7_R5_O5_RS7.txt", "Data_test_N7_R5_O6_RS7.txt", "instance_N100_R100_O100_RS25.txt",
    "instance_N100_R100_O150_RS25.txt", "instance_N100_R50_O50_RS25.txt", "instance_N200_R100_O100_RS25.txt", "instance_N200_R50_O50_RS25.txt"]

file = 4

function main(file)
    println("Reading da ta from ", data_files[file])
    data = readData("../data/" * data_files[file])
    println("Data read successfully.")


    # Set values for P, Capa, FO, and SO
    nb_picker = 5
    nb_FO = 5
    capa = 12

    data.P = nb_picker
    data.Capa = [capa for _ in 1:data.P]
    data.FO = [i for i in 1:nb_FO]
    data.SO = [i for i in 1:data.O if i ∉ data.FO]

    start_time = time()

    #opt, val, varx, vary = plneP1(data.N, data.R, data.O, data.RS, data.P, data.S, data.Q, data.Capa, data.FO, data.SO, 3, 600)
    opt, val, varx, vary, obj, LB = dantzigWolfe(data.N, data.R, data.O, data.RS, data.P, data.S, data.Q, data.Capa, data.FO, data.SO, 60)

    if opt
        println()
        println("Optimal solution found.")
        println("Objective value: ", val)
        println("x: ", varx)
        println("y: ", vary)
        println("nb of x", sum(varx))
        println("nb of y", sum(vary))
        println("time: ", time() - start_time, " seconds.")
        println("Number of iterations: ", length(obj))
        println("Lower bounds: ", LB)
        println("Objective values: ", obj)
    else
        println("No optimal solution found within the time limit.")
        println("time: ", time() - start_time, " seconds.")
    end
end

