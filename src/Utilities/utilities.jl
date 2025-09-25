"""
```julia
@ivy(x::Expr)
```

Return the expression `x` with `@inbounds` and `@views` in front of it.
"""
macro ivy end

macro ivy(x::Expr)
    return esc(quote
        @inbounds @views $x
    end)
end


# Used for testing
examples_dir() = pkgdir(PinCFlow,"examples/submit")
export examples_dir

function compute_norms(state)
	(; rho, u, v, w, pip, rhop) = state.variables.predictands
	(; jac) = state.grid
    
	l2_norm = [sqrt(sum(abs2.(rho) .* jac));
		    sqrt(sum(abs2.(u) .* jac));
    		    sqrt(sum(abs2.(v) .* jac));
		    sqrt(sum(abs2.(w) .* jac));
		    sqrt(sum(abs2.(pip) .* jac));
		sqrt(sum(abs2.(rhop) .* jac))] 
   
	linf_norm = [ maximum(abs.(rho));
                 	maximum(abs.(u));
			maximum(abs.(v));
                 	maximum(abs.(w));
                 	maximum(abs.(pip));
	maximum(abs.(rhop))];
	return l2_norm, linf_norm	
end

function analysis_callback(state)
	l2_norm, linf_norm = compute_norms(state)
	return (; l2 = l2_norm, linf = linf_norm)
end
