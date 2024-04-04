using JuMP
using CPLEX

function dantzigWolfe(N, R, O, RS, P, S, Q, Capa, FO, SO, time_limit=60)
    start_time = time()

    obj = []
    LB = []

    M = 1000000
    epsilon = 0.00000000001
    iter = 0
    updated = true
    current_x = []
    current_y = []
    current_u = []
    current_v = []

    println(Q)
    println(S)

    while updated && time() - start_time < time_limit
        updated = false

        # JuMP model with the CPLEX optimizer, number of threads, and time limit set
        model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 4, "CPX_PARAM_TILIM" => time_limit))
        set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0) # Remove the solver output

        # Decision variables
        @variable(model, lambda1[k=1:iter] >= 0)
        @variable(model, lambda2[k=1:iter] >= 0)
        @variable(model, xa >= 0)

        # Objective function
        @objective(model, Min, M * xa + (length(SO) + 1) * sum(lambda2[k] * sum(current_u[k][r] for r = 1:R) for k = 1:iter) - sum(lambda1[k] * sum(current_v[k][o] for o in SO) for k = 1:iter))

        # Constraints
        @constraint(model, con[p=1:P, i=1:N], sum(lambda1[k] * sum(current_x[k][p, o] * Q[i][o] for o = 1:O) for k = 1:iter) <= sum(lambda2[k] * sum(current_y[k][p, r] * S[i][r] for r = 1:R) for k = 1:iter))
        @constraint(model, con1, sum(lambda1[k] for k = 1:iter) + xa == 1)
        @constraint(model, con2, sum(lambda2[k] for k = 1:iter) + xa == 1)

        optimize!(model)
        #println(model)
        println("objective value = ", objective_value(model))
        push!(obj, objective_value(model))

        alpha = -JuMP.dual.(con)
        n1 = JuMP.dual.(con1)
        n2 = JuMP.dual.(con2)
        #println("alpha = ", alpha)
        #println("n1 = ", n1)
        #println("n2 = ", n2)

        subP1 = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 1, "CPX_PARAM_TILIM" => time_limit))
        set_optimizer_attribute(subP1, "CPX_PARAM_SCRIND", 0) # Remove the solver output

        # Decision variables
        @variable(subP1, x[p=1:P, o=1:O], Bin)
        @variable(subP1, 0 <= v[o=1:O] <= 1)

        # Objective function
        @objective(subP1, Min, -sum(v[o] for o in SO) + sum(sum(x[p, o] * sum(alpha[p, i] * Q[i][o] for i = 1:N) for p = 1:P) for o in 1:O))

        # Constraints
        @constraint(subP1, [o in SO], sum(x[p, o] for p = 1:P) == v[o])
        @constraint(subP1, [o in FO], sum(x[p, o] for p = 1:P) == 1)
        @constraint(subP1, [p = 1:P], sum(x[p, o] for o = 1:O) <= Capa[p])

        optimize!(subP1)
        #println("sp1 x = ", value.(x))
        #println("sp1 v = ", value.(v))
        #println("obj val of sp1 = ", objective_value(subP1), "  n1 = ", n1)
        if termination_status(subP1) == MOI.OPTIMAL
            if objective_value(subP1) - n1 <= -epsilon
                new_x = value.(x)
                new_v = value.(v)
                updated = true
                push!(current_x, new_x)
                push!(current_v, new_v)
            else
                new_x = zeros(P, O)
                for i in FO
                    for p in 1:P
                        new_x[p, i] = 1
                    end
                end
                new_v = zeros(O)
                push!(current_x, new_x)
                push!(current_v, new_v)
            end
        end


        subP2 = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 1, "CPX_PARAM_TILIM" => time_limit))
        set_optimizer_attribute(subP2, "CPX_PARAM_SCRIND", 0) # Remove the solver output

        # Decision variables
        @variable(subP2, 0 <= u[r=1:R] <= 1)
        @variable(subP2, y[p=1:P, r=1:R], Bin)

        # Objective function
        @objective(subP2, Min, (length(SO) + 1) * sum(u[r] for r = 1:R) + sum(sum(y[p, r] * sum(S[i][r] * (-alpha[p, i]) for i = 1:N) for p = 1:P) for r = 1:R))

        # Constraints
        @constraint(subP2, [r = 1:R], sum(y[p, r] for p = 1:P) == u[r])

        optimize!(subP2)
        #println("sp2 y = ", value.(y))
        #println("sp2 u = ", value.(u))
        #println("obj val of sp2 = ", objective_value(subP2), "  n2 = ", n2)
        if termination_status(subP2) == MOI.OPTIMAL
            if objective_value(subP2) - n2 <= -epsilon
                new_y = value.(y)
                new_u = value.(u)
                updated = true
                push!(current_y, new_y)
                push!(current_u, new_u)
            else
                new_y = zeros(P, R)
                new_u = zeros(R)
                push!(current_y, new_y)
                push!(current_u, new_u)
            end
        end

        println("sp1 val = ", objective_value(subP1), "  sp2 val = ", objective_value(subP2), "   sum = ", objective_value(subP1) + objective_value(subP2))
        push!(LB, objective_value(subP1) + objective_value(subP2))
        if updated
            iter += 1
        end
    end

    # JuMP model with the CPLEX optimizer, number of threads, and time limit set
    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 4, "CPX_PARAM_TILIM" => time_limit))
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0) # Remove the solver output

    # Decision variables
    @variable(model, lambda1[k=1:iter], Bin)
    @variable(model, lambda2[k=1:iter], Bin)
    @variable(model, xa >= 0)

    # Objective function
    @objective(model, Min, M * xa + (length(SO) + 1) * sum(lambda2[k] * sum(current_u[k][r] for r = 1:R) for k = 1:iter) - sum(lambda1[k] * sum(current_v[k][o] for o in SO) for k = 1:iter))

    # Constraints
    @constraint(model, con[p=1:P, i=1:N], sum(lambda1[k] * sum(current_x[k][p, o] * Q[i][o] for o = 1:O) for k = 1:iter) <= sum(lambda2[k] * sum(current_y[k][p, r] * S[i][r] for r = 1:R) for k = 1:iter))
    @constraint(model, con1, sum(lambda1[k] for k = 1:iter) + xa == 1)
    @constraint(model, con2, sum(lambda2[k] for k = 1:iter) + xa == 1)

    optimize!(model)

    # Results
    if termination_status(model) == MOI.OPTIMAL
        x = zeros(P, O)
        y = zeros(P, R)
        lambda1 = value.(lambda1)
        lambda2 = value.(lambda2)
        for i in 1:iter
            if lambda1[i] > 0
                println("x = ", current_x[i])
                x += current_x[i] * lambda1[i]
                println("x = ", x)
            end
        end
        for i in 1:iter
            if lambda2[i] > 0
                y += current_y[i] * lambda2[i]
                println("y = ", y)
            end
        end
        println("iter = ", iter)
        for i in 1:iter
            println("x_", i, " = ", current_x[i])
            println("y_", i, " = ", current_y[i])
        end
        return true, objective_value(model), x, y, obj, LB
    else
        #println("No optimal solution found within the time limit.")
        if has_values(model)
            #println("Best solution found:")
            #println("Objective value: ", objective_value(model))
            # Access the values of decision variables if needed
            # x_optimal = value.(x)
            # y_optimal = value.(y)
            return false, objective_value(model), value.(lambda1), value.(lambda2), obj, LB
        else
            return false, Inf, [], []
        end
    end
end