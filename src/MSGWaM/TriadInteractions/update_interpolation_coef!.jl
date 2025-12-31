function update_interpolation_coef end

function update_interpolation_coef!(spec_tend::TriadTendencies, nk)
    
    (; kp, m, kpl, ml, lambdakp, lambdam) = spec_tend.spec_grid

    (; c_o, alphakp, alpham, beta) = spec_tend.interp_coef


   


    for ih = 2:kpl, iz = 2:ml
        n1 = nk[ih-1, iz-1]
        n2 = nk[ih, iz-1]
        n3 = nk[ih, iz]
        n4 = nk[ih-1, iz]
        b = (-n1 + n2 - n3 + n4) / ((lambdakp - 1.0) * (lambdam - 1.0))

        beta[ih, iz] = b / (kp[ih-1] * m[iz-1])
        alphakp[ih, iz] = ((n1 - n2) / (lambdakp - 1.0) - b) / kp[ih-1]
        alpham[ih, iz] = ((n1 - n4) / (lambdam - 1.0) - b) / kp[iz-1]
        c_o[ih, iz] = n1 + (n1 - n2) / (lambdakp - 1.0) + (n1 - n4) / (lambdam - 1.0) - b
    end

    for iz = 2:ml
        n1 = max(0.0, (kp[2] * nk[1,iz-1] - kp[1] * nk[2,iz-1]) / (kp[2]-kp[1]) )
        n2 = nk[1,iz-1]
        n3 = nk[1,iz]
        n4 = max(0.0, (kp[2] * nk[1,iz] - kp[1] * nk[2,iz]) / (kp[2]-kp[1]) )
        b = (-n1 + n2 - n3 + n4) / (kp[1] * (m[iz] - m[iz-1]))
        
        beta[1, iz] = b
        alphakp[1, iz] = (n1 - n2) / kp[1] - b * m[iz-1] 
        alpham[1, iz] = (n1 - n4) / (m[iz] - m[iz-1])
        c_o[1, iz] = n1 + alpham[1, iz] * m[iz-1]
    end
    
    for ih = 2:kpl
        n1 = max(0.0, (m[2] * nk[ih-1,1] - m[1] * nk[ih-1,2]) / (m[2]-m[1]) )
        n2 = max(0.0, (m[2] * nk[ih,1] - m[1] * nk[ih,2]) / (m[2]-m[1]) )
        n3 = nk[ih,1]
        n4 = nk[ih-1,1]
        b = (-n1 + n2 - n3 + n4) / (m[1] * (kp[ih] - kp[ih-1]))

        beta[ih, 1] = b
        alphakp[ih, 1] = (n1 - n2) / (kp[ih] - kp[ih-1])
        alpham[ih, 1] = (n1 - n4) / m[1] - b * kp[ih-1] 
        c_o[ih, 1] = n1 + alphakp[ih, 1] * kp[ih-1]
    end

    # For ih = iz = 1 (linear nk = c0 - αh kh - αz kz if nks ≥ 0)
    n2 = max(0.0, (m[2] * nk[1,1] - m[1] * nk[1,2]) / (m[2]-m[1]) )
    n3 = nk[1,1]
    n4 = max(0.0, (kp[2] * nk[1,1] - kp[1] * nk[2,1]) / (kp[2]-kp[1]) )
    n1 = max(0.0, n2 - n3 + n4) 
    b = (-n1 + n2 - n3 + n4) / (m[1] * kp[1]) 

    beta[1, 1] = b
    alphakp[1, 1] = (n1 - n2) / kp[1]
    alpham[1, 1] = (n1 - n4) / m[1]
    c_o[1, 1] = n1
end