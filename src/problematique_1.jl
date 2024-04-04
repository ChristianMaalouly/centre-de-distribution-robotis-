using JuMP
using CPLEX

function plneP1(N, R, O, RS, P, S, Q, Capa, FO, SO, obj=3, time_limit=60)

    M = 1000

    # JuMP model with the CPLEX optimizer, number of threads, and time limit set
    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 8, "CPX_PARAM_TILIM" => time_limit))
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0) # Remove the solver output

    # Decision variables
    @variable(model, y[p=1:P, r=1:R], Bin)
    @variable(model, x[p=1:P, o=1:O], Bin)
    @variable(model, 0 <= u[r=1:R] <= 1)
    @variable(model, 0 <= v[o=1:O] <= 1)

    # Objective function
    if obj == 1
        @objective(model, Min, sum(u[r] for r = 1:R))
    elseif obj == 2
        @objective(model, Max, (sum(v[o] for o = 1:O)))
    elseif obj == 3
        @objective(model, Min, sum(u[r] * (length(SO) + 1) for r = 1:R) - sum(v[o] for o in SO))
    elseif obj == 4
        # error with abs
        @objective(model, Min, sum(sum((sum(x[p1, r] for r in 1:R) - sum(x[p2, r] for r in 1:R))^2 for p2 = p1+1:P) for p1 = 1:P-1))
    end

    # Constraints
    @constraint(model, [o in SO], sum(x[p, o] for p = 1:P) == v[o])
    @constraint(model, [o in FO], sum(x[p, o] for p = 1:P) == 1)
    @constraint(model, [r = 1:R], sum(y[p, r] for p = 1:P) == u[r])
    @constraint(model, [p = 1:P], sum(x[p, o] for o = 1:O) <= Capa[p])
    @constraint(model, [p = 1:P, i = 1:N], sum(x[p, o] * Q[i][o] for o = 1:O) <= sum(y[p, r] * S[i][r] for r = 1:R))


    # Solve the optimization problem
    optimize!(model)
    # Results
    if termination_status(model) == MOI.OPTIMAL
        #println("Optimal solution found.")
        #println("Objective value: ", objective_value(model))
        # Access the values of decision variables if needed
        # x_optimal = value.(x)
        # y_optimal = value.(y)
        return true, objective_value(model), value.(x), value.(y)
    else
        #println("No optimal solution found within the time limit.")
        if has_values(model)
            #println("Best solution found:")
            #println("Objective value: ", objective_value(model))
            # Access the values of decision variables if needed
            # x_optimal = value.(x)
            # y_optimal = value.(y)
            return false, objective_value(model), value.(x), value.(y)
        else
            return false, Inf, [], []
        end
    end
end