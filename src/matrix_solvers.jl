using ConjugateGradientsGPU

struct CGSolver{RealT <: Real} <: AbstractMatrixSolver
    maxiter::Int
    tol::RealT
end

function CGSolver(; maxiter = 100, tol = 1.0f-6)
    CGSolver(maxiter, tol)
end
